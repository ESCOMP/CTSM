module surfrdMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in surface data file and determining
  ! subgrid weights
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod    , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use clm_varpar      , only : nlevsoifl
  use landunit_varcon , only : numurbl
  use clm_varcon      , only : grlnd
  use clm_varctl      , only : iulog
  use clm_varctl      , only : use_cndv, use_crop, use_fates
  use surfrdUtilsMod  , only : check_sums_equal_1, apply_convert_ocean_to_land, collapse_crop_types
  use surfrdUtilsMod  , only : collapse_to_dominant, collapse_crop_var, collapse_individual_lunits
  use ncdio_pio       , only : file_desc_t, var_desc_t, ncd_pio_openfile, ncd_pio_closefile
  use ncdio_pio       , only : ncd_io, check_var, ncd_inqfdims, check_dim_size, ncd_inqdid, ncd_inqdlen
  use pio
  use spmdMod
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: surfrd_compat_check     ! Check that this surface dataset is compatible
  public :: surfrd_get_data         ! Read surface dataset and determine subgrid weights
  public :: surfrd_get_num_patches  ! Read surface dataset to determine maxsoil_patches and numcft
  public :: surfrd_get_nlevurb      ! Read surface dataset to determine nlevurb

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: surfrd_special   ! Read the special landunits
  private :: surfrd_veg_all   ! Read all of the vegetated landunits
  private :: surfrd_veg_dgvm  ! Read vegetated landunits for DGVM mode
  private :: surfrd_pftformat ! Read crop pfts in file format where they are part of the vegetated land unit
  private :: surfrd_cftformat ! Read crop pfts in file format where they are on their own landunit
  !
  ! !PRIVATE DATA MEMBERS:
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  subroutine check_domain_attributes(ncid, begg, endg, ldomain, info)
    ! !DESCRIPTION:
    ! Checks for mismatches between the land domain and a surface or similar dataset's domain.
    !
    ! !USES:
    use domainMod, only : domain_type, domain_init, domain_clean
    !
    ! !ARGUMENTS
    type(file_desc_t), intent(inout) :: ncid ! netcdf id for input file
    integer, intent(in) :: begg, endg
    type(domain_type), intent(in) :: ldomain  ! land domain
    character(len=*), intent(in) :: info  ! information to include in messages
    !
    ! !LOCAL VARIABLES
    type(domain_type) :: inputdata_domain  ! local domain associated with input dataset
    logical :: readvar  ! true => variable is on dataset
    logical :: istype_domain  ! true => input file is of type domain
    character(len=16) :: lon_var, lat_var  ! names of lat/lon on dataset
    logical :: isgrid2d  ! true => input grid is 2d
    integer :: ni, nj, ns  ! domain sizes
    integer :: n
    real(r8) :: rmaxlon, rmaxlat  ! local min/max vars

   character(len=32) :: subname = 'check_domain_attributes'  ! subroutine name

    call check_var(ncid=ncid, varname='xc', readvar=readvar)
    if (readvar) then
       istype_domain = .true.
    else
       call check_var(ncid=ncid, varname='LONGXY', readvar=readvar)
       if (readvar) then
          istype_domain = .false.
       else
          call endrun( msg=' ERROR: unknown '//info//' domain type---'//errMsg(sourcefile, __LINE__))
       end if
    end if
    if (istype_domain) then
       lon_var  = 'xc'
       lat_var  = 'yc'
    else
       lon_var  = 'LONGXY'
       lat_var  = 'LATIXY'
    end if
    if ( masterproc )then
       write(iulog,*) trim(subname),' ',info,' lon_var = ',trim(lon_var),' lat_var =',trim(lat_var)
    end if

    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    call domain_init(inputdata_domain, isgrid2d, ni, nj, begg, endg, subgrid_level=grlnd)

    call ncd_io(ncid=ncid, varname=lon_var, flag='read', data=inputdata_domain%lonc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: lon var NOT on '//info//' dataset---'//errMsg(sourcefile, __LINE__))

    call ncd_io(ncid=ncid, varname=lat_var, flag='read', data=inputdata_domain%latc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: lat var NOT on '//info//' dataset---'//errMsg(sourcefile, __LINE__))

    rmaxlon = 0.0_r8
    rmaxlat = 0.0_r8
    do n = begg,endg
       if (ldomain%lonc(n)-inputdata_domain%lonc(n) > 300.) then
          rmaxlon = max(rmaxlon,abs(ldomain%lonc(n)-inputdata_domain%lonc(n)-360._r8))
       elseif (ldomain%lonc(n)-inputdata_domain%lonc(n) < -300.) then
          rmaxlon = max(rmaxlon,abs(ldomain%lonc(n)-inputdata_domain%lonc(n)+360._r8))
       else
          rmaxlon = max(rmaxlon,abs(ldomain%lonc(n)-inputdata_domain%lonc(n)))
       endif
       rmaxlat = max(rmaxlat,abs(ldomain%latc(n)-inputdata_domain%latc(n)))
    enddo
    if (rmaxlon > 0.001_r8 .or. rmaxlat > 0.001_r8) then
       write(iulog,*)' ERROR: '//info//' dataset vs. land domain lon/lat mismatch error', rmaxlon,rmaxlat
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    call domain_clean(inputdata_domain)
  end subroutine check_domain_attributes

  !-----------------------------------------------------------------------
  subroutine surfrd_compat_check ( lfsurdat )
    !
    ! !DESCRIPTION:
    ! Check compatability for this surface dataset and abort with an error if it's not
    !
    ! !USES:
    use ncdio_pio, only : check_att
    ! !ARGUMENTS:
    character(len=*), intent(in) :: lfsurdat    ! surface dataset filename
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid  ! netcdf id
    logical :: exists          ! If attribute or variable was found on the file
    integer :: status          ! Status return code
    real(r4) :: version        ! Version number on the dataset
    ! NOTE: Only increment the expected_version when surface datasets are incompatible with the previous version
    !       If datasets are just updated data and backwards compatble leave the expected version alone
    real(r4), parameter :: expected_version = 5.3_r4
    character(len=50) :: description
    character(len=*), parameter :: version_name = 'Dataset_Version'

    call ncd_pio_openfile (ncid, trim(lfsurdat), 0)
    call check_att(ncid, pio_global, version_name, exists)
    if (exists) then
       status = pio_get_att(ncid, pio_global, version_name, version)
    else
       ! For a few previous versions guess on the compatability version based on existence of variables
       call check_var( ncid, 'PCT_OCEAN', exists)
       if (exists) then
         version = 5.2_r4
       else
         call check_var( ncid, 'CONST_HARVEST_SH1', exists)
         if (exists) then
            version = 5.0_r4
         else
            call check_var( ncid, 'GLACIER_REGION', exists)
            if (exists) then
               version = 4.5_r4
            else
               ! This is a version before the main clm4_5 dataseta so marking it as 0 for unknown
               version = 0.0_r4
            end if
         end if
       end if
    end if
    call ncd_pio_closefile(ncid)
    if ( (version /= expected_version) )then
      if ( version < expected_version )then
         description = 'older'
         if ( version == 0.0_r4 ) description = trim(description)//' than 4.5'
      else if ( version > expected_version )then
         description = 'newer'
      end if
      if ( masterproc )then
         write(iulog,*) 'Input surface dataset is: ', trim(lfsurdat)
         write(iulog,'(3a)') 'This surface dataset is ', trim(description), ' and incompatible with this version of CTSM'
         write(iulog,'(a,f3.1,a,f3.1)') 'Dataset version = ', version, ' Version expected = ', expected_version
         write(iulog,*) errMsg(sourcefile, __LINE__)
      end if
      call endrun(msg="ERROR: Incompatible surface dataset")
    end if

  end subroutine surfrd_compat_check

  !-----------------------------------------------------------------------
  subroutine surfrd_get_data (begg, endg, ldomain, lfsurdat, lhillslope_file, actual_numcft)
    !
    ! !DESCRIPTION:
    ! Read the surface dataset and create subgrid weights.
    ! The model's surface dataset recognizes 6 basic land cover types within a grid
    ! cell: lake, wetland, urban, glacier, glacier_mec and vegetated. The vegetated
    ! part of the grid cell cosists of up to [maxsoil_patches] patches. These
    ! subgrid patches are read in explicitly for each grid cell. This is in
    ! contrast to LSMv1, where the patches were built implicitly from biome types.
    !    o real latitude  of grid cell (degrees)
    !    o real longitude of grid cell (degrees)
    !    o integer surface type: 0 = ocean or 1 = land
    !    o integer soil color (1 to 20) for use with soil albedos
    !    o real soil texture, %sand, for thermal and hydraulic properties
    !    o real soil texture, %clay, for thermal and hydraulic properties
    !    o real % of cell covered by lake    for use as subgrid patch
    !    o real % of cell covered by wetland for use as subgrid patch
    !    o real % of cell that is urban      for use as subgrid patch
    !    o real % of cell that is glacier    for use as subgrid patch
    !    o real % of cell that is glacier_mec for use as subgrid patch
    !    o integer PFTs
    !    o real % abundance PFTs (as a percent of vegetated area)
    !
    ! !USES:
    use clm_varctl  , only : create_crop_landunit, convert_ocean_to_land, collapse_urban, &
                             toosmall_soil, toosmall_crop, toosmall_glacier, &
                             toosmall_lake, toosmall_wetland, toosmall_urban, &
                             n_dom_landunits, &
                             use_hillslope
    use fileutils           , only : getfil
    use domainMod           , only : domain_type
    use clm_instur          , only : wt_lunit, topo_glc_mec, pct_urban_max
    use landunit_varcon     , only : max_lunit, istsoil, isturb_MIN, isturb_MAX
    use dynSubgridControlMod, only : get_flanduse_timeseries
    use dynSubgridControlMod, only : get_do_transient_lakes
    use dynSubgridControlMod, only : get_do_transient_urban

    !
    ! !ARGUMENTS:
    integer,          intent(in) :: begg, endg, actual_numcft
    type(domain_type),intent(in) :: ldomain     ! land domain
    character(len=*), intent(in) :: lfsurdat    ! surface dataset filename
    character(len=*), intent(in) :: lhillslope_file ! hillslope dataset filename
    !
    ! !LOCAL VARIABLES:
    character(len=256):: locfn                ! local file name
    integer, parameter :: n_dom_urban = 1     ! # of dominant urban landunits
    type(file_desc_t) :: ncid                 ! netcdf id for lfsurdat
    type(file_desc_t) :: ncid_hillslope       ! netcdf id for lhillslope_file

    character(len=32) :: subname = 'surfrd_get_data'    ! subroutine name
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read surface boundary data .....'
       if (lfsurdat == ' ') then
          write(iulog,*)'lfsurdat must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
       if (use_hillslope .and. lhillslope_file == ' ') then
          write(iulog,*)'lhillslope_file must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    endif

    wt_lunit(:,:) = 0._r8
    topo_glc_mec(:,:) = 0._r8

    ! Read surface data

    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    if (use_hillslope) then
       call getfil( lhillslope_file, locfn, 0 )
       call ncd_pio_openfile (ncid_hillslope, trim(locfn), 0)
    end if

    ! Compare dataset domain attributes to ldomain attributes
    call check_domain_attributes(ncid, begg, endg, ldomain, 'surface')
    if (use_hillslope) then
       call check_domain_attributes(ncid_hillslope, begg, endg, ldomain, 'hillslope')
    end if

    ! Obtain special landunit info

    call surfrd_special(begg, endg, ncid, ldomain%ns)

    ! Obtain vegetated landunit info

    call surfrd_veg_all(begg, endg, ncid, ncid_hillslope, ldomain%ns, actual_numcft)

    if (use_cndv) then
       call surfrd_veg_dgvm(begg, endg)
    end if

    call ncd_pio_closefile(ncid)

    call check_sums_equal_1(wt_lunit, begg, 'wt_lunit', subname)

    if (convert_ocean_to_land) then
       call apply_convert_ocean_to_land(wt_lunit(begg:endg,:), begg, endg)
    end if

    ! if collapse_urban = .true.
    ! collapse urban landunits to the dominant urban landunit

    if (collapse_urban) then
       call collapse_to_dominant(wt_lunit(begg:endg,isturb_MIN:isturb_MAX), isturb_MIN, isturb_MAX, begg, endg, n_dom_urban)
    end if

    ! Select N dominant landunits
    ! ---------------------------
    ! n_dom_landunits set by user in namelist
    ! Call resembles the surfrd_veg_all call to the same subr that selects
    ! n_dom_pfts (also set by the user in the namelist)

    call collapse_to_dominant(wt_lunit(begg:endg,:), istsoil, max_lunit, &
                              begg, endg, n_dom_landunits)

    ! Remove landunits using thresholds set by user in namelist
    ! ---------------------------------------------------------
    ! Thresholds are set in the namelist parameters toosmall_* in units of %.
    ! TODO Remove corresponding thresholds from the mksurfdat tool
    !      Found 2 such cases (had expected to encounter one per landunit):
    !         mkurbanparCommonMod.F90 MIN_DENS = 0.1 and
    !         mksurfdat.F90 toosmallPFT = 1.e-10

    call collapse_individual_lunits(wt_lunit, begg, endg, toosmall_soil, &
                                    toosmall_crop, toosmall_glacier, &
                                    toosmall_lake, toosmall_wetland, &
                                    toosmall_urban)

    if ( masterproc )then
       write(iulog,*) 'Successfully read surface boundary data'
       write(iulog,*)
    end if

    ! read the lakemask (necessary for initialisation of dynamical lakes)
    if (get_do_transient_lakes()) then
        call surfrd_lakemask(begg, endg)
    end if

    ! read the urbanmask (necessary for initialization of dynamical urban)
    if (get_do_transient_urban()) then
        call surfrd_urbanmask(begg, endg)
    else
        ! Set this to zero here. pct_urban_max is used in subgridWeightsMod to check 
        ! whether urban landunits should be run virtually.
        pct_urban_max(:,:) = 0._r8
    end if
    
  end subroutine surfrd_get_data

!-----------------------------------------------------------------------

  subroutine surfrd_get_num_patches (lfsurdat, actual_maxsoil_patches, actual_numpft, actual_numcft)
    !
    ! !DESCRIPTION:
    ! Read maxsoil_patches and numcft from the surface dataset
    !
    ! !USES:
    use fileutils   , only : getfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: lfsurdat  ! surface dataset filename
    integer, intent(out) :: actual_maxsoil_patches  ! value from surface dataset
    integer, intent(out) :: actual_numcft           ! cft value from sfc dataset
    integer, intent(out) :: actual_numpft           ! pft value from sfc dataset
    
    !
    ! !LOCAL VARIABLES:
    character(len=256):: locfn                ! local file name
    type(file_desc_t) :: ncid                 ! netcdf file id
    integer :: dimid                          ! netCDF dimension id
    logical :: cft_dim_exists                 ! dimension exists on dataset
    integer :: check_numpft                   ! Surface dataset count of numpft, should
                                              ! match maxsoil_patches - actual_numcft
    character(len=32) :: subname = 'surfrd_get_num_patches'  ! subroutine name
    !-----------------------------------------------------------------------

    if (masterproc) then
       if(use_fates)then
          write(iulog,*) 'Attempting to read numcft from the surface data .....'
       else
          write(iulog,*) 'Attempting to read maxsoil_patches and numcft from the surface data .....'
       end if
       if (lfsurdat == ' ') then
          write(iulog,*)'lfsurdat must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    endif

    ! Open surface dataset
    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Read numcft
    call ncd_inqdid(ncid, 'cft', dimid, cft_dim_exists)
    if ( cft_dim_exists ) then
       call ncd_inqdlen(ncid, dimid, actual_numcft, 'cft')
    else
       actual_numcft = 0
    end if
    
    ! Read maxsoil_patches
    call ncd_inqdlen(ncid, dimid, actual_maxsoil_patches, 'lsmpft')
    actual_numpft = actual_maxsoil_patches - actual_numcft

    call ncd_inqdlen(ncid, dimid, check_numpft, 'natpft')

    if(check_numpft.ne.actual_numpft)then
       write(iulog,*)'the sum of the cftdim and the natpft dim should match the lsmpft dim in the surface file'
       write(iulog,*)'natpft: ',check_numpft
       write(iulog,*)'lsmpft: ',actual_maxsoil_patches
       write(iulog,*)'cft: ',actual_numcft
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    
    if ( masterproc )then
       write(iulog,*) 'Successfully read maxsoil_patches and numcft from the surface data'
       write(iulog,*)
    end if

  end subroutine surfrd_get_num_patches

!-----------------------------------------------------------------------
  subroutine surfrd_get_nlevurb (lfsurdat, actual_nlevurb)
    !
    ! !DESCRIPTION:
    ! Read nlevurb from the surface dataset
    !
    ! !USES:
    use fileutils   , only : getfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: lfsurdat  ! surface dataset filename
    integer, intent(out) :: actual_nlevurb    ! nlevurb from surface dataset
    !
    ! !LOCAL VARIABLES:
    character(len=256):: locfn                ! local file name
    type(file_desc_t) :: ncid                 ! netcdf file id
    integer :: dimid                          ! netCDF dimension id
    character(len=32) :: subname = 'surfrd_get_nlevurb'  ! subroutine name
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read nlevurb from the surface data .....'
       if (lfsurdat == ' ') then
          write(iulog,*)'lfsurdat must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    endif

    ! Open surface dataset
    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Read nlevurb
    call ncd_inqdlen(ncid, dimid, actual_nlevurb, 'nlevurb')

    if ( masterproc )then
       write(iulog,*) 'Successfully read nlevurb from the surface data'
       write(iulog,*)
    end if

  end subroutine surfrd_get_nlevurb

!-----------------------------------------------------------------------
  subroutine surfrd_special(begg, endg, ncid, ns)
    !
    ! !DESCRIPTION:
    ! Determine weight with respect to gridcell of all special "patches" as well
    ! as soil color and percent sand and clay
    !
    ! !USES:
    use clm_varpar      , only : maxpatch_glc, nlevurb
    use landunit_varcon , only : isturb_MIN, isturb_MAX, istdlak, istwet, istice, istocn
    use clm_instur      , only : wt_lunit, urban_valid, wt_glc_mec, topo_glc_mec
    use UrbanParamsType , only : CheckUrban
    !
    ! !ARGUMENTS:
    integer          , intent(in)    :: begg, endg
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    integer          , intent(in)    :: ns     ! domain size
    !
    ! !LOCAL VARIABLES:
    integer          :: n,nl,nurb,g   ! indices
    integer          :: dimid,varid   ! netCDF id's
    real(r8)         :: nlevsoidata(nlevsoifl)
    logical          :: found         ! temporary for error check
    integer          :: nindx         ! temporary for error check
    integer          :: ier           ! error status
    logical          :: readvar
    real(r8),pointer :: pctgla(:)     ! percent of grid cell is glacier
    real(r8),pointer :: pctlak(:)     ! percent of grid cell is lake
    real(r8),pointer :: pctwet(:)     ! percent of grid cell is wetland
    real(r8),pointer :: pctocn(:)     ! percent of grid cell is ocean
    real(r8),pointer :: pcturb(:,:)   ! percent of grid cell is urbanized
    integer ,pointer :: urban_region_id(:)
    real(r8),pointer :: pcturb_tot(:) ! percent of grid cell is urban (sum over density classes)
    real(r8),pointer :: pctspec(:)    ! percent of spec lunits wrt gcell
    integer          :: dens_index    ! urban density index
    real(r8)         :: closelat,closelon
    integer, parameter :: urban_invalid_region = 0   ! urban_region_id indicating invalid point
    character(len=32) :: subname = 'surfrd_special'  ! subroutine name
!-----------------------------------------------------------------------

    allocate(pctgla(begg:endg))
    allocate(pctlak(begg:endg))
    allocate(pctwet(begg:endg))
    allocate(pctocn(begg:endg))
    allocate(pcturb(begg:endg,numurbl))
    allocate(pcturb_tot(begg:endg))
    allocate(urban_region_id(begg:endg))
    allocate(pctspec(begg:endg))

    call check_dim_size(ncid, 'nlevsoi', nlevsoifl)

       ! Obtain non-grid surface properties of surface dataset other than percent patch

    call ncd_io(ncid=ncid, varname='PCT_OCEAN', flag='read', data=pctocn, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg= &
       ' ERROR: PCT_OCEAN NOT on surfdata file but required when running ctsm5.2 or newer; ' // &
       ' you are advised to generate a new surfdata file using the mksurfdata_esmf tool ' // errMsg(sourcefile, __LINE__))

    call ncd_io(ncid=ncid, varname='PCT_WETLAND', flag='read', data=pctwet, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_WETLAND  NOT on surfdata file'//errMsg(sourcefile, __LINE__))

    call ncd_io(ncid=ncid, varname='PCT_LAKE'   , flag='read', data=pctlak, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_LAKE NOT on surfdata file'//errMsg(sourcefile, __LINE__))

    call ncd_io(ncid=ncid, varname='PCT_GLACIER', flag='read', data=pctgla, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_GLACIER NOT on surfdata file'//errMsg(sourcefile, __LINE__))

    ! Read urban info
    if (nlevurb == 0) then
      ! If PCT_URBAN is not multi-density then set pcturb to zero
      pcturb = 0._r8
      urban_valid(begg:endg) = .false.
      write(iulog,*)'PCT_URBAN is not multi-density, pcturb set to 0'
    else
      call ncd_io(ncid=ncid, varname='PCT_URBAN'  , flag='read', data=pcturb, &
           dim1name=grlnd, readvar=readvar)
      if (.not. readvar) call endrun( msg=' ERROR: PCT_URBAN NOT on surfdata file'//errMsg(sourcefile, __LINE__))

      call ncd_io(ncid=ncid, varname='URBAN_REGION_ID', flag='read', data=urban_region_id, &
           dim1name=grlnd, readvar=readvar)
      if (.not. readvar) call endrun( msg= ' ERROR: URBAN_REGION_ID NOT on surfdata file'//errMsg(sourcefile, __LINE__))
      where (urban_region_id == urban_invalid_region)
         urban_valid = .false.
      elsewhere
         urban_valid = .true.
      end where
    end if
    if ( nlevurb == 0 )then
       if ( any(pcturb > 0.0_r8) ) then
          call endrun( msg=' ERROR: PCT_URBAN MUST be zero when nlevurb=0'//errMsg(sourcefile, __LINE__))
       end if
    end if

    pcturb_tot(:) = 0._r8
    do n = 1, numurbl
       do nl = begg,endg
          pcturb_tot(nl) = pcturb_tot(nl) + pcturb(nl,n)
       enddo
    enddo

    ! Read glacier info

    call check_dim_size(ncid, 'nglcec',   maxpatch_glc   )
    call check_dim_size(ncid, 'nglcecp1', maxpatch_glc+1 )

    call ncd_io(ncid=ncid, varname='PCT_GLC_MEC', flag='read', data=wt_glc_mec, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_GLC_MEC NOT on surfdata file'//errMsg(sourcefile, __LINE__))

    wt_glc_mec(:,:) = wt_glc_mec(:,:) / 100._r8
    call check_sums_equal_1(wt_glc_mec, begg, 'wt_glc_mec', subname)

    call ncd_io(ncid=ncid, varname='TOPO_GLC_MEC',  flag='read', data=topo_glc_mec, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: TOPO_GLC_MEC NOT on surfdata file'//errMsg(sourcefile, __LINE__))

    topo_glc_mec(:,:) = max(topo_glc_mec(:,:), 0._r8)

    pctspec = pctwet + pctlak + pcturb_tot + pctgla + pctocn

    ! Error check: sum of glacier, lake, wetland, urban, ocean must be < 100

    found = .false.
    do nl = begg,endg
       if (pctspec(nl) > 100._r8+1.e-04_r8) then
          found = .true.
          nindx = nl
          exit
       end if
       if (found) exit
    end do
    if ( found ) then
       write(iulog,*)'surfrd error: patch cover>100 for nl=',nindx
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Determine wt_lunit for special landunits

    do nl = begg,endg

       wt_lunit(nl,istdlak) = pctlak(nl) / 100._r8

       ! Until ctsm5.1 we would label ocean points as wetland in fsurdat
       ! files. Starting with ctsm5.2 we label ocean points as ocean
       ! (always 100%) and wetland points as wetland.
       wt_lunit(nl,istwet) = pctwet(nl) / 100._r8
       wt_lunit(nl,istocn) = pctocn(nl) / 100._r8
       wt_lunit(nl,istice) = pctgla(nl) / 100._r8

       do n = isturb_MIN, isturb_MAX
          dens_index = n - isturb_MIN + 1
          wt_lunit(nl,n) = pcturb(nl,dens_index) / 100._r8
       end do

    end do

    call CheckUrban(begg, endg, pcturb(begg:endg,:), subname)

    deallocate(pctgla,pctlak,pctwet,pctocn,pcturb,pcturb_tot,urban_region_id,pctspec)

  end subroutine surfrd_special

  !-----------------------------------------------------------------------
  subroutine surfrd_cftformat( ncid, begg, endg, wt_cft, fert_cft, cftsize, natpft_size )
    !
    ! !DESCRIPTION:
    !     Handle generic crop types for file format where they are on their own
    !     crop landunit and read in as Crop Function Types.
    ! !USES:
    use clm_instur      , only : wt_nat_patch, irrig_method
    use clm_varpar      , only : cft_size, cft_lb, natpft_lb, cft_ub, natpft_ub
    use IrrigationMod   , only : irrig_method_unset
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid         ! netcdf id
    integer          , intent(in)    :: begg, endg
    integer          , intent(in)    :: cftsize      ! CFT size
    real(r8), pointer, intent(inout) :: wt_cft(:,:)  ! CFT weights
    real(r8), pointer, intent(inout) :: fert_cft(:,:)! Fertilizer
    integer          , intent(in)    :: natpft_size  ! natural PFT size
    !
    ! !LOCAL VARIABLES:
    logical  :: readvar                              ! is variable on dataset
    real(r8),pointer :: array2D(:,:)                 ! local array
    character(len=32) :: subname = 'surfrd_cftformat'! subroutine name
!-----------------------------------------------------------------------
    SHR_ASSERT_ALL_FL((lbound(wt_cft)          == (/begg, cft_lb/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wt_cft, dim=1)   == (/endg/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wt_cft, dim=2)   >= (/cftsize+1-cft_lb/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((lbound(fert_cft)        == (/begg, cft_lb/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fert_cft, dim=1) == (/endg/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fert_cft, dim=2) >= (/cftsize+1-cft_lb/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wt_nat_patch)    >= (/endg,natpft_size-1+natpft_lb/)), sourcefile, __LINE__)

    call check_dim_size(ncid, 'cft',    cftsize)
    call check_dim_size(ncid, 'natpft', natpft_size)
    
    call ncd_io(ncid=ncid, varname='PCT_CFT', flag='read', data=wt_cft, &
            dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_CFT NOT on surfdata file'//errMsg(sourcefile, __LINE__))

    if ( cft_size > 0 )then
       call ncd_io(ncid=ncid, varname='CONST_FERTNITRO_CFT', flag='read', data=fert_cft, &
               dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if ( masterproc ) &
                write(iulog,*) ' WARNING: CONST_FERTNITRO_CFT NOT on surfdata file zero out'
          fert_cft = 0.0_r8
       end if
    else
       fert_cft = 0.0_r8
    end if

    if ( cft_size > 0 )then
       call ncd_io(ncid=ncid, varname='irrigation_method', flag='read', data=irrig_method, &
               dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if ( masterproc ) &
                write(iulog,*) ' WARNING: irrigation_method NOT on surfdata file; using default'
          irrig_method = irrig_method_unset
       end if
    else
       irrig_method = irrig_method_unset
    end if

    allocate( array2D(begg:endg,1:natpft_size) )
    call ncd_io(ncid=ncid, varname='PCT_NAT_PFT', flag='read', data=array2D, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_NAT_PFT NOT on surfdata file'//errMsg(sourcefile, __LINE__))
    wt_nat_patch(begg:,natpft_lb:natpft_size-1+natpft_lb) = array2D(begg:,:)
    deallocate( array2D )
 

  end subroutine surfrd_cftformat

  subroutine surfrd_wtfates( ncid, begg, endg )

    !--------------------------------------------------------------------------
    !     This routine evaluates the natural and crop functional
    !     type fractions in the surface file and returns them to
    !     a single, concatenated vector.  These weights
    !     are only used for a satellite phenology run.
    !     Note that FATES will actually allocate a different number of patches
    !     and will use a mapping table to connect its own pft and cft
    !     definitions to those it finds in the surface file.
    !--------------------------------------------------------------------------
    
    ! !USES:
    use clm_instur      , only : wt_nat_patch, wt_lunit
    use clm_varpar      , only : cft_size, surfpft_lb, surfpft_ub
    use landunit_varcon , only : istsoil, istcrop
    
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid         ! netcdf id
    integer          , intent(in)    :: begg, endg

    !
    ! !LOCAL VARIABLES:
    logical  :: readvar                              ! is variable on dataset
    real(r8),pointer :: array2d_pft(:,:)                 ! local array
    real(r8),pointer :: array2d_cft(:,:)                 ! local array
    integer :: g,p
    integer :: cft_dimlen,natpft_dimlen,dimid
    
    character(len=32) :: subname = 'surfrd_fates'! subroutine name
    
    call ncd_inqdlen(ncid, dimid, cft_dimlen, 'cft')
    call ncd_inqdlen(ncid, dimid, natpft_dimlen, 'natpft')

    ! double check that cft_dimlen+natpft_dimlen = natpft_size
    if((cft_dimlen+natpft_dimlen).ne.(surfpft_ub-surfpft_lb+1))then
       call endrun( msg=' ERROR: PCT+CFT dimlen does not match array size for wt_nat_patch when fates is on'//errMsg(sourcefile, __LINE__))
    end if
    
    allocate( array2d_cft(begg:endg,1:cft_dimlen) )
    allocate( array2d_pft(begg:endg,1:natpft_dimlen) )
    
    call ncd_io(ncid=ncid, varname='PCT_CFT', flag='read', data=array2d_cft, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_CFT NOT on surfdata file'//errMsg(sourcefile, __LINE__))
    
    call ncd_io(ncid=ncid, varname='PCT_NAT_PFT', flag='read', data=array2d_pft, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_NAT_PFT NOT on surfdata file'//errMsg(sourcefile, __LINE__))

    ! In fates, all the weights in both the cft and pfts go into this array
    ! It is only used by SP mode, and it can choose what PFTs to align with
    
    wt_nat_patch(begg:,0:natpft_dimlen-1) = array2d_pft(begg:,:)
    wt_nat_patch(begg:,natpft_dimlen:natpft_dimlen+cft_dimlen-1) = array2d_cft(begg:,:)

    if(.false.)then
       ! Scale the weights by the lu weights from the dataset
       do g = begg, endg
          
          if((wt_lunit(g,istsoil)+wt_lunit(g,istcrop))>0._r8)then
             do p = 0,natpft_dimlen-1
                wt_nat_patch(g,p) = wt_nat_patch(g,p) * wt_lunit(g,istsoil)/(wt_lunit(g,istsoil)+wt_lunit(g,istcrop))
             end do
             do p = natpft_dimlen,natpft_dimlen+cft_dimlen-1
                wt_nat_patch(g,p) = wt_nat_patch(g,p) * wt_lunit(g,istcrop)/(wt_lunit(g,istsoil)+wt_lunit(g,istcrop))
             end do
          else
             wt_nat_patch(g,:) = 0._r8
             wt_nat_patch(g,0) = 100._r8
          end if
       end do
       ! Add the crop weight to the natveg weight and zero the crop weight
       wt_lunit(begg:,istsoil) = wt_lunit(begg:,istsoil) + wt_lunit(begg:,istcrop)
       wt_lunit(begg:,istcrop) = 0._r8
       
    else
       ! Legacy method
       do g = begg, endg
          if ( wt_lunit(g,istcrop) > 0.0_r8 )then
             ! Move CFT over to PFT and do weighted average of the crop and soil parts
             wt_nat_patch(g,0:natpft_dimlen-1) = wt_nat_patch(g,0:natpft_dimlen-1) * wt_lunit(g,istsoil)
             wt_nat_patch(g,natpft_dimlen:natpft_dimlen+cft_dimlen-1)       = &
                  wt_nat_patch(g,natpft_dimlen:natpft_dimlen+cft_dimlen-1) * wt_lunit(g,istcrop)
             wt_lunit(g,istsoil) = (wt_lunit(g,istsoil) + wt_lunit(g,istcrop)) ! Add crop landunit to soil landunit
             wt_nat_patch(g,:)   =  wt_nat_patch(g,:) / wt_lunit(g,istsoil)
             wt_lunit(g,istcrop) = 0.0_r8                ! Zero out crop CFT's
          else
             wt_nat_patch(g,natpft_dimlen:natpft_dimlen+cft_dimlen-1) = 0.0_r8    ! Make sure generic crops are zeroed out
          end if
       end do
       
    end if
    
    deallocate(array2d_cft,array2d_pft)
    

  end subroutine surfrd_wtfates

  
  !-----------------------------------------------------------------------
  
  subroutine surfrd_pftformat( begg, endg, ncid )
    !
    ! !DESCRIPTION:
    !     Handle generic crop types for file format where they are part of the
    !     natural vegetation landunit.
    ! !USES:
    use clm_instur      , only : fert_cft, irrig_method, wt_nat_patch
    use clm_varpar      , only : natpft_size, cft_size, natpft_lb
    use IrrigationMod   , only : irrig_method_unset
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begg, endg
    type(file_desc_t), intent(inout) :: ncid                    ! netcdf id
    !
    ! !LOCAL VARIABLES:
    logical  :: cft_dim_exists                 ! does the dimension 'cft' exist on the dataset?
    integer  :: dimid                          ! netCDF id's
    logical  :: readvar                        ! is variable on dataset
    character(len=32) :: subname = 'surfrd_pftformat'! subroutine name
!-----------------------------------------------------------------------
    SHR_ASSERT_ALL_FL((ubound(wt_nat_patch) == (/endg, natpft_size-1+natpft_lb/)), sourcefile, __LINE__)

    call check_dim_size(ncid, 'natpft', natpft_size)
    
    ! If cft_size == 0, then we expect to be running with a surface dataset
    ! that does
    ! NOT have a PCT_CFT array (or CONST_FERTNITRO_CFT array), and thus does not have a 'cft' dimension.
    ! Make sure
    ! that's the case.
    call ncd_inqdid(ncid, 'cft', dimid, cft_dim_exists)
    if (cft_dim_exists) then
       call endrun( msg= ' ERROR: unexpectedly found cft dimension on dataset when cft_size=0'// &
               ' (if the surface dataset has a separate crop landunit, then the code'// &
               ' must also have a separate crop landunit, and vice versa)'//&
               errMsg(sourcefile, __LINE__))
    end if
    call ncd_io(ncid=ncid, varname='CONST_FERTNITRO_CFT', flag='read', data=fert_cft, &
            dim1name=grlnd, readvar=readvar)
    if (readvar) then
       call endrun( msg= ' ERROR: unexpectedly found CONST_FERTNITRO_CFT on dataset when cft_size=0'// &
               ' (if the surface dataset has a separate crop landunit, then the code'// &
               ' must also have a separate crop landunit, and vice versa)'//&
               errMsg(sourcefile, __LINE__))
    end if
    fert_cft = 0.0_r8

    call ncd_io(ncid=ncid, varname='irrigation_method', flag='read', data=irrig_method, &
            dim1name=grlnd, readvar=readvar)
    if (readvar) then
       call endrun( msg= ' ERROR: unexpectedly found irrigation_method on dataset when cft_size=0'// &
               ' (if the surface dataset has a separate crop landunit, then the code'// &
               ' must also have a separate crop landunit, and vice versa)'//&
               errMsg(sourcefile, __LINE__))
    end if
    irrig_method = irrig_method_unset

    call ncd_io(ncid=ncid, varname='PCT_NAT_PFT', flag='read', data=wt_nat_patch, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_NAT_PFT NOT on surfdata file'//errMsg(sourcefile, __LINE__))
       
  end subroutine surfrd_pftformat

!-----------------------------------------------------------------------
  subroutine surfrd_veg_all(begg, endg, ncid, ncid_hillslope, ns, actual_numcft)
    !
    ! !DESCRIPTION:
    ! Determine weight arrays for non-dynamic landuse mode
    !
    ! !USES:
    use clm_varctl      , only : create_crop_landunit, use_fates, n_dom_pfts, use_hillslope
    use clm_varpar      , only : natpft_lb, natpft_ub, natpft_size, cft_size, cft_lb, cft_ub
    use clm_varpar      , only : surfpft_lb, surfpft_ub
    use clm_instur      , only : wt_lunit, wt_nat_patch, wt_cft, fert_cft
    use landunit_varcon , only : istsoil, istcrop
    use surfrdUtilsMod  , only : convert_cft_to_pft

    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begg, endg, actual_numcft
    type(file_desc_t),intent(inout) :: ncid   ! netcdf id for fsurdat
    type(file_desc_t),intent(inout) :: ncid_hillslope ! netcdf id for hillslope_file
    integer          ,intent(in)    :: ns     ! domain size
    !
    ! !LOCAL VARIABLES:
    integer  :: dimid                                ! netCDF id's
    integer  :: cftsize                              ! size of CFT's
    logical  :: readvar                              ! is variable on dataset
    logical  :: cft_dim_exists                       ! does the dimension 'cft' exist on the dataset?
    real(r8),pointer :: arrayl(:)                    ! local array
    real(r8),pointer :: array2DCFT(:,:)              ! local 2D array for CFTs
    real(r8),pointer :: array2DFERT(:,:)             ! local 2D array for fertilizer
    character(len=32) :: subname = 'surfrd_veg_all'  ! subroutine name

    !-----------------------------------------------------------------------
    !
    ! Read in variables that are handled the same for all formats
    !

    ! This temporary array is needed because ncd_io expects a pointer, so we can't
    ! directly pass wt_lunit(begg:endg,istsoil)
    allocate(arrayl(begg:endg))

    call ncd_io(ncid=ncid, varname='PCT_NATVEG', flag='read', data=arrayl, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_NATVEG NOT on surfdata file'//errMsg(sourcefile, __LINE__))
    wt_lunit(begg:endg,istsoil) = arrayl(begg:endg)

    call ncd_io(ncid=ncid, varname='PCT_CROP', flag='read', data=arrayl, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_CROP NOT on surfdata file'//errMsg(sourcefile, __LINE__))
    wt_lunit(begg:endg,istcrop) = arrayl(begg:endg)

    deallocate(arrayl)

    
    ! Check the file format for CFT's and handle accordingly
    if ( actual_numcft > 0 ) then
       if ( create_crop_landunit )then
          ! Cases: Full crop in file and in model
          !        Generic crop in file and in model
          !        Full crop in file, generic crop in model
          cftsize = actual_numcft
          allocate(array2DCFT (begg:endg,cft_lb:cftsize-1+cft_lb))
          allocate(array2DFERT(begg:endg,cft_lb:cftsize-1+cft_lb))
          call surfrd_cftformat( ncid, begg, endg, array2DCFT, array2DFERT, cftsize, natpft_size )
          wt_cft  (begg:,cft_lb:) = array2DCFT (begg:,cft_lb:cft_ub)
          fert_cft(begg:,cft_lb:) = array2DFERT(begg:,cft_lb:cft_ub)
          deallocate(array2DCFT)
          deallocate(array2DFERT)
       else if ( .not. create_crop_landunit )then
          if ( masterproc ) write(iulog,*) "WARNING: New CFT-based format surface datasets should be run with ", &
                                           "create_crop_landunit=T"
          if ( use_fates ) then
             if ( masterproc ) write(iulog,*) "WARNING: When fates is on we allow new CFT based surface datasets ", &
                                              "to be used with create_crop_land FALSE"

             call surfrd_wtfates( ncid, begg, endg )
             ! Set the weighting on the crop patches to zero
             fert_cft(begg:,cft_lb:cft_ub) = 0.0_r8
             wt_cft(begg:,cft_lb:cft_ub) = 0.0_r8

          else
             call endrun( msg=' ERROR: New format surface datasets require create_crop_landunit TRUE'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun( msg=' ERROR: Problem figuring out how to handle new format input fsurdat file'//errMsg(sourcefile, __LINE__))
       end if
    else if ( (.not. cft_dim_exists) .and. (.not. create_crop_landunit) )then
       if ( masterproc ) write(iulog,*) "WARNING: The PFT format is an unsupported format that will be removed in the future!"
       ! Check dimension size
       call surfrd_pftformat( begg, endg, ncid )                                 ! Format where crop is part of the natural veg. landunit
    else
       call endrun( msg=' ERROR: Problem figuring out format of input fsurdat file'//errMsg(sourcefile, __LINE__))
    end if

    ! Do some checking

    if ( (cft_size == 0) .and. any(wt_lunit(begg:endg,istcrop) > 0._r8) ) then
       call endrun( msg=' ERROR: if PCT_CROP > 0 anywhere, then cft_size must be > 0'// &
               ' (if the surface dataset has a separate crop landunit, then the code'// &
               ' must also have a separate crop landunit, and vice versa)'//&
               errMsg(sourcefile, __LINE__))
    end if

    ! Obtain hillslope hydrology information and modify pft weights
    if (use_hillslope) then
       call surfrd_hillslope(begg, endg, ncid_hillslope, ns)
    endif

    ! Convert from percent to fraction
    wt_lunit(begg:endg,istsoil) = wt_lunit(begg:endg,istsoil) / 100._r8
    wt_lunit(begg:endg,istcrop) = wt_lunit(begg:endg,istcrop) / 100._r8
    wt_nat_patch(begg:endg,:)   = wt_nat_patch(begg:endg,:) / 100._r8
    wt_cft(begg:endg,:)         = wt_cft(begg:endg,:) / 100._r8

    ! Check sum of vegetation adds to 1
    call check_sums_equal_1(wt_nat_patch, begg, 'wt_nat_patch', subname)

    ! if ( use_fates ) wt_cft = 0 because called convert_cft_to_pft, else...
    if ( .not. use_fates ) then
       ! Check sum of vegetation adds to 1
       call check_sums_equal_1(wt_cft, begg, 'wt_cft', subname)
    end if

    ! Call collapse_crop_types: allows need to maintain only 78-pft input data
    ! For use_crop = .false. collapsing 78->16 pfts or 16->16 or some new
       !    configuration
    ! For use_crop = .true. most likely collapsing 78 to the list of crops for
    !    which the CLM includes parameterizations
    ! The call collapse_crop_types also appears in subroutine dyncrop_interp
    call collapse_crop_types(wt_cft(begg:endg,:), fert_cft(begg:endg,:), cft_size, begg, endg, verbose=.true.)
    
    ! Collapse crop variables as needed
    ! The call to collapse_crop_var also appears in subroutine dyncrop_interp
    ! - fert_cft TODO Is this call redundant because it simply sets the crop
    !                 variable to 0 where is_pft_known_to_model = .false.?
    call collapse_crop_var(fert_cft(begg:endg,:), cft_size, begg, endg)

    ! Call collapse_to_dominant: enhance ctsm performance with fewer active pfts
    ! Collapsing to the top N dominant pfts (n_dom_pfts set in namelist).
    ! - Bare ground could be up to 1 patch before collapsing.
    ! - Pfts could be up to 14 before collapsing if create_crop_landunit = .T.
    ! - Pfts could be up to 16 before collapsing if create_crop_landunit = .F.
    ! TODO Add the same call to subroutine dynpft_interp for transient runs
    
    call collapse_to_dominant(wt_nat_patch(begg:endg,:), surfpft_lb, surfpft_ub, &
         begg, endg, n_dom_pfts)
    
  end subroutine surfrd_veg_all

  !-----------------------------------------------------------------------
  subroutine surfrd_veg_dgvm(begg, endg)
    !
    ! !DESCRIPTION:
    ! Determine weights for CNDV mode.
    !
    ! !USES:
    use pftconMod , only : noveg
    use clm_instur, only : wt_nat_patch
    !
    ! !ARGUMENTS:
    integer, intent(in) :: begg, endg
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'surfrd_veg_dgvm'
    !-----------------------------------------------------------------------

    ! Bare ground gets 100% weight; all other natural patches are zeroed out
    wt_nat_patch(begg:endg, :)     = 0._r8
    wt_nat_patch(begg:endg, noveg) = 1._r8

    call check_sums_equal_1(wt_nat_patch, begg, 'wt_nat_patch', subname)

  end subroutine surfrd_veg_dgvm

  !-----------------------------------------------------------------------
  subroutine surfrd_hillslope(begg, endg, ncid, ns)
    !
    ! !DESCRIPTION:
    ! Determine number of hillslopes and columns for hillslope hydrology mode
    !
    ! !USES:
    use clm_instur, only : ncolumns_hillslope, wt_nat_patch
    use clm_varctl, only : nhillslope,max_columns_hillslope
    use clm_varpar, only : natpft_size, natpft_lb, natpft_ub
    use ncdio_pio,  only : ncd_inqdid, ncd_inqdlen
    use pftconMod , only : noveg
    use HillslopeHydrologyMod, only : pft_distribution_method, pft_standard, pft_from_file, pft_uniform_dominant_pft, pft_lowland_dominant_pft, pft_lowland_upland
    use array_utils, only: find_k_max_indices
    use surfrdUtilsMod, only: collapse_to_dominant

    !
    ! !ARGUMENTS:
    integer, intent(in) :: begg, endg
    type(file_desc_t),intent(inout) :: ncid   ! netcdf id
    integer          ,intent(in)    :: ns     ! domain size
    !
    ! !LOCAL VARIABLES:
    integer  :: g, nh, m, n                       ! index
    integer  :: dimid,varid                    ! netCDF id's
    integer  :: ier                            ! error status
    integer, allocatable  :: max_indices(:)    ! largest weight pft indices
    logical  :: readvar                        ! is variable on dataset
    integer,pointer :: arrayl(:)               ! local array (needed because ncd_io expects a pointer)
    character(len=32) :: subname = 'surfrd_hillslope'  ! subroutine name
    logical, allocatable :: do_not_collapse(:)
    integer :: n_dominant
    !-----------------------------------------------------------------------

    ! number of hillslopes per landunit
    call ncd_inqdid(ncid,'nhillslope',dimid,readvar)
    if (.not. readvar) then
       call endrun( msg=' ERROR: nhillslope not on surface data file'//errMsg(sourcefile, __LINE__))
    else
       call ncd_inqdlen(ncid,dimid,nh)
       nhillslope = nh
    endif
    ! maximum number of columns per landunit
    call ncd_inqdid(ncid,'nmaxhillcol',dimid,readvar)
    if (.not. readvar) then
       call endrun( msg=' ERROR: nmaxhillcol not on surface data file'//errMsg(sourcefile, __LINE__))
    else
       call ncd_inqdlen(ncid,dimid,nh)
       max_columns_hillslope = nh
    endif
    ! actual number of columns per landunit
    allocate(arrayl(begg:endg))
    call ncd_io(ncid=ncid, varname='nhillcolumns', flag='read', data=arrayl, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun( msg=' ERROR: nhillcolumns not on surface data file'//errMsg(sourcefile, __LINE__))
    else
       ncolumns_hillslope(begg:endg) = arrayl(begg:endg)
    endif
    deallocate(arrayl)

    ! pft_from_file and pft_lowland_upland assume that 1 pft
    ! will exist on each hillslope column.  In prepration, set one
    ! pft weight to 100 and the rest to 0.  The vegetation type
    ! (patch%itype) will be reassigned when initHillslope is called later.
    if(pft_distribution_method == pft_from_file .or. &
         pft_distribution_method == pft_lowland_upland) then
       do g = begg, endg
          ! If hillslopes will be used in a gridcell, modify wt_nat_patch, otherwise use original patch distribution
          if(ncolumns_hillslope(g) > 0) then
             ! First patch gets 100% weight; all other natural patches are zeroed out
             wt_nat_patch(g,:)         = 0._r8
             wt_nat_patch(g,natpft_lb) = 100._r8
          endif
       enddo

    else if (pft_distribution_method == pft_uniform_dominant_pft &
        .or. pft_distribution_method == pft_lowland_dominant_pft) then

       ! If hillslopes will be used in a gridcell, modify wt_nat_patch,
       ! otherwise use original patch distribution
       allocate(do_not_collapse(begg:endg))
       do_not_collapse(begg:endg) = .false.
       do g = begg, endg
          if (ncolumns_hillslope(g) == 0) then
             do_not_collapse(g) = .true.
          end if
       end do

       if (pft_distribution_method == pft_uniform_dominant_pft) then
         ! pft_uniform_dominant_pft uses the patch with the
         ! largest weight for all hillslope columns in the gridcell
         n_dominant = 1
       else if (pft_distribution_method == pft_lowland_dominant_pft) then
         ! pft_lowland_dominant_pft uses the two patches with the
         ! largest weights for the hillslope columns in the gridcell
         n_dominant = 2
       else
          call endrun( msg=' ERROR: unrecognized hillslope_pft_distribution_method'//errMsg(sourcefile, __LINE__))
       end if

       call collapse_to_dominant(wt_nat_patch(begg:endg,:), natpft_lb, natpft_ub, begg, endg, n_dominant, do_not_collapse)
       deallocate(do_not_collapse)

    else if (pft_distribution_method /= pft_standard) then
      call endrun( msg=' ERROR: unrecognized hillslope_pft_distribution_method'//errMsg(sourcefile, __LINE__))
    endif

  end subroutine surfrd_hillslope

  subroutine surfrd_lakemask(begg, endg)
    !
    ! !DESCRIPTION:
    ! Reads the lake mask, indicating where lakes are and will grow
    ! of the landuse.timeseries file.
    ! Necessary for the initialisation of the lake land units
    !
    ! !USES:
     use clm_instur           , only : pct_lake_max
     use dynSubgridControlMod , only : get_flanduse_timeseries
     use clm_varctl           , only : fname_len
     use fileutils            , only : getfil
    !
    ! !ARGUMENTS:
    integer,           intent(in)    :: begg, endg
    !
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t)         :: ncid_dynuse          ! netcdf id for landuse timeseries file
    character(len=256)        :: locfn                ! local file name
    character(len=fname_len)  :: fdynuse              ! landuse.timeseries filename
    logical                   :: readvar
    !
    character(len=*), parameter :: subname = 'surfrd_lakemask'
    !
    !-----------------------------------------------------------------------

    ! get filename of landuse_timeseries file
    fdynuse = get_flanduse_timeseries()

    if (masterproc) then
       write(iulog,*) 'Attempting to read landuse.timeseries data .....'
       if (fdynuse == ' ') then
          write(iulog,*)'fdynuse must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    call getfil(fdynuse, locfn, 0 )

   ! open landuse_timeseries file
    call ncd_pio_openfile (ncid_dynuse, trim(locfn), 0)

    ! read the lakemask
    call ncd_io(ncid=ncid_dynuse, varname='PCT_LAKE_MAX'  , flag='read', data=pct_lake_max, &
           dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_LAKE_MAX is not on landuse.timeseries file'//errMsg(sourcefile, __LINE__))

    ! close landuse_timeseries file again
    call ncd_pio_closefile(ncid_dynuse)

  end subroutine surfrd_lakemask

  !-----------------------------------------------------------------------
  subroutine surfrd_urbanmask(begg, endg)
    !
    ! !DESCRIPTION:
    ! Reads the urban mask, indicating where urban areas are and will grow
    ! of the landuse.timeseries file.
    ! Necessary for the initialization of the urban land units.
    ! All urban density types will intialize if any type exists or will grow.
    !
    ! !USES:
     use clm_instur           , only : pct_urban_max
     use dynSubgridControlMod , only : get_flanduse_timeseries
     use clm_varctl           , only : fname_len
     use fileutils            , only : getfil
    !
    ! !ARGUMENTS:
    integer,           intent(in)    :: begg, endg
    !
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t)         :: ncid_dynuse          ! netcdf id for landuse timeseries file
    character(len=256)        :: locfn                ! local file name
    character(len=fname_len)  :: fdynuse              ! landuse.timeseries filename
    logical                   :: readvar
    !
    character(len=*), parameter :: subname = 'surfrd_urbanmask'
    !
    !-----------------------------------------------------------------------

    ! get filename of landuse_timeseries file
    fdynuse = get_flanduse_timeseries()

    if (masterproc) then
       write(iulog,*) 'Attempting to read landuse.timeseries data .....'
       if (fdynuse == ' ') then
          write(iulog,*)'fdynuse must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    call getfil(fdynuse, locfn, 0 )

   ! open landuse_timeseries file
    call ncd_pio_openfile (ncid_dynuse, trim(locfn), 0)

    ! read the urbanmask
    call ncd_io(ncid=ncid_dynuse, varname='PCT_URBAN_MAX', flag='read', data=pct_urban_max, &
           dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_URBAN_MAX is not on landuse.timeseries file'//errMsg(sourcefile, __LINE__))

    ! close landuse_timeseries file again
    call ncd_pio_closefile(ncid_dynuse)

  end subroutine surfrd_urbanmask
  
end module surfrdMod
