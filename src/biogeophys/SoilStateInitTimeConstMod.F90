module SoilStateInitTimeConstMod

  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Set hydraulic and thermal properties 
  !
  ! !USES
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use SoilStateType , only : soilstate_type
  use LandunitType  , only : lun                
  use ColumnType    , only : col                
  use PatchType     , only : patch                
  use abortUtils    , only : endrun
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: SoilStateInitTimeConst
  public  :: readParams

  ! PRIVATE FUNCTIONS MADE PUBLIC Just for unit-testing:
  public :: ThresholdSoilMoistZender2003
  public :: ThresholdSoilMoistKok2014
  public :: MassFracClay
  public :: MassFracClayLeung2023
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ReadNL
  !
  ! !PUBLIC DATA:
  real(r8), public :: organic_max  ! organic matter (kg/m3) where soil is assumed to act like peat

  ! !PRIVATE DATA:
  type, private :: params_type
     real(r8) :: tkd_sand            ! Thermal conductivity of sand (W/m/K)
     real(r8) :: tkd_clay            ! Thermal conductivity of clay (W/m/K)
     real(r8) :: tkd_om              ! Thermal conductivity of dry organic matter (Farouki, 1981) (W/m/K)
     real(r8) :: tkm_om              ! Thermal conductivity of organic matter (Farouki, 1986) (W/m/K)
     real(r8) :: pd                  ! Particle density of soil (kg/m3)
     real(r8) :: csol_clay           ! Heat capacity of clay *10^6 (J/K/m3)
     real(r8) :: csol_om             ! Heat capacity of peat soil *10^6 (Farouki, 1986) (J/K/m3)
     real(r8) :: csol_sand           ! Heat capacity of sand *10^6 (J/K/m3)
     real(r8) :: bsw_sf              ! Scale factor for bsw (unitless)
     real(r8) :: hksat_sf            ! Scale factor for hksat (unitless)
     real(r8) :: sucsat_sf           ! Scale factor for sucsat (unitless)
     real(r8) :: watsat_sf           ! Scale factor for watsat (unitless)
     real(r8) :: sand_pf             ! Perturbation factor (via addition) for percent sand (percent)
     real(r8) :: clay_pf             ! Perturbation factor (via addition) for percent clay of clay+silt (percent)
     real(r8) :: om_frac_sf          ! Scale factor for organic matter fraction (unitless)
  end type params_type
  type(params_type), private ::  params_inst

  ! Control variables (from namelist)
  logical, private :: organic_frac_squared ! If organic fraction should be squared (as in CLM4.5)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------
  !
contains

  !-----------------------------------------------------------------------
  subroutine ReadNL( nlfilename )
    !
    ! !DESCRIPTION:
    ! Read namelist for SoilStateType
    !
    ! !USES:
    use shr_mpi_mod    , only : shr_mpi_bcast
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use fileutils      , only : getavu, relavu, opnfil
    use clm_nlUtilsMod , only : find_nlgroup_name
    use clm_varctl     , only : iulog
    use spmdMod        , only : mpicom, masterproc
    use abortUtils     , only : endrun    
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: nlfilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'SoilState_readnl'  ! subroutine name
    !-----------------------------------------------------------------------

    character(len=*), parameter :: nl_name  = 'clm_soilstate_inparm'  ! Namelist name
                                                                      ! MUST agree with name in namelist and read
    namelist / clm_soilstate_inparm / organic_frac_squared

    ! preset values

    organic_frac_squared = .false.

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in '//nl_name//' namelist'
       call opnfil (nlfilename, unitn, 'F')
       call find_nlgroup_name(unitn, nl_name, status=ierr)
       if (ierr == 0) then
          read(unit=unitn, nml=clm_soilstate_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading '//nl_name//' namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding '//nl_name//' namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )

    end if

    call shr_mpi_bcast(organic_frac_squared, mpicom)

  end subroutine ReadNL

  !-----------------------------------------------------------------------
  subroutine readParams( ncid )
    !
    ! !USES:
    use ncdio_pio, only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'readParams_SoilStateInitTimeConst'
    !--------------------------------------------------------------------

    ! Thermal conductivity of sand (W/m/K)
    call readNcdioScalar(ncid, 'tkd_sand', subname, params_inst%tkd_sand)
    ! Thermal conductivity of clay (W/m/K)
    call readNcdioScalar(ncid, 'tkd_clay', subname, params_inst%tkd_clay)
    ! Thermal conductivity of dry organic matter (Farouki, 1981) (W/m/K)
    call readNcdioScalar(ncid, 'tkd_om', subname, params_inst%tkd_om)
    ! Thermal conductivity of organic matter (Farouki, 1986) (W/m/K)
    call readNcdioScalar(ncid, 'tkm_om', subname, params_inst%tkm_om)
    ! Particle density of soil (kg/m3)
    call readNcdioScalar(ncid, 'pd', subname, params_inst%pd)
    ! Heat capacity of clay *10^6 (J/K/m3)
    call readNcdioScalar(ncid, 'csol_clay', subname, params_inst%csol_clay)
    ! Heat capacity of peat soil *10^6 (Farouki, 1986) (J/K/m3)
    call readNcdioScalar(ncid, 'csol_om', subname, params_inst%csol_om)
    ! Heat capacity of sand *10^6 (J/K/m3)
    call readNcdioScalar(ncid, 'csol_sand', subname, params_inst%csol_sand)
    ! Scale factor for bsw (unitless)
    call readNcdioScalar(ncid, 'bsw_sf', subname, params_inst%bsw_sf)
    ! Scale factor for hksat (unitless)
    call readNcdioScalar(ncid, 'hksat_sf', subname, params_inst%hksat_sf)
    ! Scale factor for sucsat (unitless)
    call readNcdioScalar(ncid, 'sucsat_sf', subname, params_inst%sucsat_sf)
    ! Scale factor for watsat (unitless)
    call readNcdioScalar(ncid, 'watsat_sf', subname, params_inst%watsat_sf)
    ! Perturbation factor (via addition) for percent sand (percent)
    call readNcdioScalar(ncid, 'sand_pf', subname, params_inst%sand_pf)
    ! Perturbation factor  (via addition) for percent clay of clay+silt (percent)
    call readNcdioScalar(ncid, 'clay_pf', subname, params_inst%clay_pf)
    ! Scale factor for organic matter fraction (unitless)
    call readNcdioScalar(ncid, 'om_frac_sf', subname, params_inst%om_frac_sf)

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilStateInitTimeConst(bounds, soilstate_inst, nlfilename) 
    !
    ! !USES:
    use shr_log_mod         , only : errMsg => shr_log_errMsg
    use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
    use decompMod           , only : bounds_type, subgrid_level_gridcell
    use abortutils          , only : endrun
    use spmdMod             , only : masterproc
    use ncdio_pio           , only : file_desc_t, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use ncdio_pio           , only : ncd_pio_openfile, ncd_pio_closefile, ncd_inqdlen
    use clm_varpar          , only : numrad
    use clm_varpar          , only : nlevsoi, nlevgrnd, nlevlak, nlevsoifl, nlayer, nlayert, nlevmaxurbgrnd, nlevsno
    use clm_varcon          , only : zsoi, dzsoi, zisoi, spval
    use clm_varcon          , only : secspday, denh2o, denice, grlnd
    use clm_varctl          , only : use_cn, use_lch4, use_fates
    use clm_varctl          , only : iulog, fsurdat, paramfile, soil_layerstruct_predefined
    use landunit_varcon     , only : istdlak, istwet, istsoil, istcrop, istice
    use column_varcon       , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv 
    use fileutils           , only : getfil
    use organicFileMod      , only : organicrd 
    use FuncPedotransferMod , only : pedotransf, get_ipedof
    use RootBiophysMod      , only : init_vegrootfr
    use GridcellType        , only : grc
    use shr_dust_emis_mod   , only : is_dust_emis_zender, is_dust_emis_leung
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds  
    type(soilstate_type) , intent(inout) :: soilstate_inst
    character(len=*)     , intent(in)    :: nlfilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer            :: p, lev, c, l, g, j            ! indices
    real(r8)           :: om_frac                       ! organic matter fraction
    real(r8)           :: om_watsat_lake = 0.9_r8       ! porosity of organic soil
    real(r8)           :: om_hksat_lake  = 0.1_r8       ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat_lake = 10.3_r8      ! saturated suction for organic matter (Letts, 2000)
    real(r8)           :: om_b_lake      = 2.7_r8       ! Clapp Hornberger paramater for oragnic soil (Letts, 2000) (lake)
    real(r8)           :: om_watsat                     ! porosity of organic soil
    real(r8)           :: om_hksat                      ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat                     ! saturated suction for organic matter (mm)(Letts, 2000)
    real(r8)           :: om_b                          ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8)           :: zsapric        = 0.5_r8       ! depth (m) that organic matter takes on characteristics of sapric peat
    real(r8)           :: pcalpha        = 0.5_r8       ! percolation threshold
    real(r8)           :: pcbeta         = 0.139_r8     ! percolation exponent
    real(r8)           :: pc_lake        = 0.5_r8       ! percolation threshold
    real(r8)           :: perc_frac                     ! "percolating" fraction of organic soil
    real(r8)           :: perc_norm                     ! normalize to 1 when 100% organic soil
    real(r8)           :: uncon_hksat                   ! series conductivity of mineral/organic soil
    real(r8)           :: uncon_frac                    ! fraction of "unconnected" soil
    real(r8)           :: bd                            ! bulk density of dry soil material [kg/m^3]
    real(r8)           :: tkm                           ! mineral conductivity
    real(r8)           :: xksat                         ! maximum hydraulic conductivity of soil [mm/s]
    real(r8)           :: clay,sand                     ! temporaries
    real(r8)           :: perturbed_sand                ! temporary for paramfile implementation of +/- sand percentage
    real(r8)           :: residual_clay_frac            ! temporary for paramfile implementation of +/- residual clay percentage
    real(r8)           :: perturbed_residual_clay_frac  ! temporary for paramfile implementation of +/- residual clay percentage
    real(r8)           :: dust_moist_fact               ! tuning factor for soil moisture effect on limiting dust emissions, used by Charlie Zender. Simone Tilmes suggested to change this parameter into a namelist variable for easier CESM tuning. dmleung added 30 Sep 2024
    integer            :: dimid                         ! dimension id
    logical            :: readvar 
    type(file_desc_t)  :: ncid                          ! netcdf id
    real(r8) ,pointer  :: zsoifl (:)                    ! Output: [real(r8) (:)]  original soil midpoint 
    real(r8) ,pointer  :: zisoifl (:)                   ! Output: [real(r8) (:)]  original soil interface depth 
    real(r8) ,pointer  :: gti (:)                       ! read in - fmax 
    real(r8) ,pointer  :: sand3d (:,:)                  ! read in - soil texture: percent sand (needs to be a pointer for use in ncdio)
    real(r8) ,pointer  :: clay3d (:,:)                  ! read in - soil texture: percent clay (needs to be a pointer for use in ncdio)
    real(r8) ,pointer  :: organic3d (:,:)               ! read in - organic matter: kg/m3 (needs to be a pointer for use in ncdio)
    character(len=256) :: locfn                         ! local filename
    integer            :: ipedof  
    integer            :: begp, endp
    integer            :: begc, endc
    integer            :: begg, endg
    integer :: found  ! flag that equals 0 if not found and 1 if found
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    do c = begc,endc
       soilstate_inst%smpmin_col(c) = -1.e8_r8
    end do

    ! --------------------------------------------------------------------
    ! Read namelist
    ! --------------------------------------------------------------------

    call ReadNL( nlfilename )

    ! --------------------------------------------------------------------
    ! Initialize root fraction (computing from surface, d is depth in meter):
    ! --------------------------------------------------------------------

    ! Currently pervious road has same properties as soil
    do c = begc,endc
       l = col%landunit(c)

       if (lun%urbpoi(l) .and. col%itype(c) == icol_road_perv) then 
          do lev = 1, nlevgrnd
             soilstate_inst%rootfr_road_perv_col(c,lev) = 0._r8
          enddo
          do lev = 1,nlevsoi
             soilstate_inst%rootfr_road_perv_col(c,lev) = 1.0_r8/real(nlevsoi,r8)
          end do
! remove roots below bedrock layer
          soilstate_inst%rootfr_road_perv_col(c,1:col%nbedrock(c)) = &
               soilstate_inst%rootfr_road_perv_col(c,1:col%nbedrock(c)) &
               + sum(soilstate_inst%rootfr_road_perv_col(c,col%nbedrock(c)+1:nlevsoi)) &
               /real(col%nbedrock(c))
          soilstate_inst%rootfr_road_perv_col(c,col%nbedrock(c)+1:nlevsoi) = 0._r8
       end if
    end do

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          soilstate_inst%rootfr_col (c,nlevsoi+1:nlevgrnd) = 0._r8
       else
          ! Inactive CH4 columns 
          ! (Also includes (lun%itype(l)==istdlak .and.  allowlakeprod), which used to be
          ! in a separate branch of the conditional)
          soilstate_inst%rootfr_col (c,:) = spval
       end if
    end do

    ! Initialize root fraction 
    ! Note that fates has its own root fraction root fraction routine and should not
    ! use the following since it depends on patch%itype - which fates should not use

    if (.not. use_fates) then
        call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
             soilstate_inst%rootfr_patch(begp:endp,1:nlevgrnd),'water')
        call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
             soilstate_inst%crootfr_patch(begp:endp,1:nlevgrnd),'carbon')
     end if

    ! --------------------------------------------------------------------
    ! dynamic memory allocation
    ! --------------------------------------------------------------------

    allocate(sand3d(begg:endg,nlevsoifl))
    allocate(clay3d(begg:endg,nlevsoifl))

    ! Determine organic_max from parameter file

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_io(ncid=ncid, varname='organic_max', flag='read', data=organic_max, readvar=readvar)
    if ( .not. readvar ) call endrun(msg=' ERROR: organic_max not on param file'//errMsg(sourcefile, __LINE__))
    call ncd_pio_closefile(ncid)

    ! --------------------------------------------------------------------
    ! Read surface dataset
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil color, sand and clay boundary data .....'
    end if

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! Read in organic matter dataset 

    allocate(organic3d(begg:endg,nlevsoifl))
    call organicrd(organic3d)

    ! Read in sand and clay data

    call ncd_io(ncid=ncid, varname='PCT_SAND', flag='read', data=sand3d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: PCT_SAND NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if

    call ncd_io(ncid=ncid, varname='PCT_CLAY', flag='read', data=clay3d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: PCT_CLAY NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if

    do p = begp,endp
       g = patch%gridcell(p)
       if ( sand3d(g,1)+clay3d(g,1) == 0.0_r8 )then
          if ( any( sand3d(g,:)+clay3d(g,:) /= 0.0_r8 ) )then
             call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, &
                  msg='found depth points that do NOT sum to zero when surface does'//&
                  errMsg(sourcefile, __LINE__)) 
          end if
          sand3d(g,:) = 1.0_r8
          clay3d(g,:) = 1.0_r8
       end if
       if ( any( sand3d(g,:)+clay3d(g,:) == 0.0_r8 ) )then
          call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, &
               msg='after setting, found points sum to zero'//errMsg(sourcefile, __LINE__))
       end if

       soilstate_inst%sandfrac_patch(p) = sand3d(g,1)/100.0_r8
       soilstate_inst%clayfrac_patch(p) = clay3d(g,1)/100.0_r8
    end do

    ! Read fmax

    allocate(gti(begg:endg))
    call ncd_io(ncid=ncid, varname='FMAX', flag='read', data=gti, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: FMAX NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if
    do c = begc, endc
       g = col%gridcell(c)
       soilstate_inst%wtfact_col(c) = gti(g)
    end do
    deallocate(gti)

    ! Close file

    call ncd_pio_closefile(ncid)

    ! --------------------------------------------------------------------
    ! get original soil depths to be used in interpolation of sand and clay
    ! --------------------------------------------------------------------

    ! Note that the depths on the file are assumed to be the same as the depths in the
    ! model when running with 10SL_3.5m. Ideally zsoifl and zisoifl would be read from
    ! the surface dataset rather than assumed here.
    !
    ! We need to specify zsoifl down to nlevsoifl+1 (rather than just nlevsoifl) so that
    ! we can get the appropriate zisoifl at level nlevsoifl (i.e., the bottom interface
    ! depth).
    allocate(zsoifl(1:nlevsoifl+1), zisoifl(0:nlevsoifl))
    do j = 1, nlevsoifl+1
       zsoifl(j) = 0.025_r8*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
    enddo

    zisoifl(0) = 0._r8
    do j = 1, nlevsoifl
       zisoifl(j) = 0.5_r8*(zsoifl(j)+zsoifl(j+1))         !interface depths
    enddo

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: non-lake
    ! --------------------------------------------------------------------

    !   urban roof, sunwall and shadewall thermal properties used to 
    !   derive thermal conductivity and heat capacity are set to special 
    !   value because thermal conductivity and heat capacity for urban 
    !   roof, sunwall and shadewall are prescribed in SoilThermProp.F90 
    !   in SoilPhysicsMod.F90

    do c = begc, endc
       g = col%gridcell(c)
       l = col%landunit(c)

       ! istwet and istice and
       ! urban roof, sunwall, shadewall properties set to special value
       if (lun%itype(l)==istwet .or. lun%itype(l)==istice .or. &
           (lun%urbpoi(l) .and. col%itype(c) /= icol_road_perv .and. &
                                col%itype(c) /= icol_road_imperv)) then

          do lev = 1,nlevmaxurbgrnd
             soilstate_inst%watsat_col(c,lev) = spval
          end do

          do lev = 1,nlevgrnd
             soilstate_inst%bsw_col(c,lev)    = spval
             soilstate_inst%watfc_col(c,lev)  = spval
             soilstate_inst%hksat_col(c,lev)  = spval
             soilstate_inst%sucsat_col(c,lev) = spval
             soilstate_inst%watdry_col(c,lev) = spval 
             soilstate_inst%watopt_col(c,lev) = spval 
             soilstate_inst%bd_col(c,lev)     = spval 
             if (lev <= nlevsoi) then
                soilstate_inst%cellsand_col(c,lev) = spval
                soilstate_inst%cellclay_col(c,lev) = spval
                soilstate_inst%cellorg_col(c,lev)  = spval
             end if
          end do

          do lev = 1,nlevgrnd
             soilstate_inst%tkmg_col(c,lev)   = spval
             soilstate_inst%tksatu_col(c,lev) = spval
             soilstate_inst%tkdry_col(c,lev)  = spval
             soilstate_inst%csol_col(c,lev)= spval
          end do

       else

          do lev = 1,nlevgrnd
             ! Top-most model soil level corresponds to dataset's top-most soil
             ! level regardless of corresponding depths
             if (lev .eq. 1) then
                clay = clay3d(g,1)
                sand = sand3d(g,1)
                om_frac = min(params_inst%om_frac_sf*organic3d(g,1)/organic_max, 1._r8)
             else if (lev <= nlevsoi) then
                found = 0  ! reset value
                if (zsoi(lev) <= zisoifl(1)) then
                   ! Search above the dataset's range of zisoifl depths
                   clay = clay3d(g,1)
                   sand = sand3d(g,1)
                   om_frac = min(params_inst%om_frac_sf*organic3d(g,1)/organic_max, 1._r8)
                   found = 1
                else if (zsoi(lev) > zisoifl(nlevsoifl)) then
                   ! Search below the dataset's range of zisoifl depths
                   clay = clay3d(g,nlevsoifl)
                   sand = sand3d(g,nlevsoifl)
                   om_frac = min(params_inst%om_frac_sf*organic3d(g,nlevsoifl)/organic_max, 1._r8)
                   found = 1
                else
                   ! For remaining model soil levels, search within dataset's
                   ! range of zisoifl values. Look for model node depths
                   ! that are between the dataset's interface depths.
                   do j = 1,nlevsoifl-1
                      if (zsoi(lev) > zisoifl(j) .AND. zsoi(lev) <= zisoifl(j+1)) then
                         clay = clay3d(g,j+1)
                         sand = sand3d(g,j+1)
                         om_frac = min(params_inst%om_frac_sf*organic3d(g,j+1)/organic_max, 1._r8)
                         found = 1
                      endif
                      if (found == 1) exit  ! no need to stay in the loop
                   end do
                end if
                ! If not found, then something's wrong
                if (found == 0) then
                   write(iulog,*) 'For model soil level =', lev
                   call endrun(msg="ERROR finding a soil dataset depth to interpolate the model depth to"//errmsg(sourcefile, __LINE__))
                end if
             else  ! if lev > nlevsoi
                clay = clay3d(g,nlevsoifl)
                sand = sand3d(g,nlevsoifl)
                om_frac = 0._r8
             endif

             if (organic_frac_squared) then
                om_frac = om_frac**2._r8
             end if

             if (lun%urbpoi(l)) then
                om_frac = 0._r8 ! No organic matter for urban
             end if

             if (lev <= nlevsoi) then
                ! This is separated into sections for non-perturbation and perturbation of sand/clay
                ! because the perturbation code is not bfb when sand_pf=clay_pf=0. This occurs because
                ! of a divide and then a multiply in the code.
                if (params_inst%sand_pf == 0._r8 .and. params_inst%clay_pf == 0._r8) then
                   soilstate_inst%cellsand_col(c,lev) = sand
                   soilstate_inst%cellclay_col(c,lev) = clay
                else
                   ! by default, will read sand and clay from the surface dataset
                   !     - sand_pf can be used to perturb the absolute percent sand
                   !     - clay_pf can be used to perturb what percent of (clay+silt) is clay
                   if (sand<100._r8) then
                      residual_clay_frac              = clay/(100._r8-sand)
                   else
                      residual_clay_frac              = 0.5_r8
                   end if
                   perturbed_sand                     = min(100._r8,max(0._r8,sand+params_inst%sand_pf))
                   perturbed_residual_clay_frac       = min(1._r8,max(0._r8,residual_clay_frac + &
                                                        params_inst%clay_pf/100._r8))
                   soilstate_inst%cellsand_col(c,lev) = perturbed_sand
                   soilstate_inst%cellclay_col(c,lev) = (100._r8-perturbed_sand)*perturbed_residual_clay_frac
                end if
                soilstate_inst%cellorg_col(c,lev)  = om_frac*organic_max
             end if

             if (lun%itype(l) /= istdlak) then  ! soil columns of both urban and non-urban types

                ! Note that the following properties are overwritten for urban impervious road 
                ! layers that are not soil in SoilThermProp.F90 within SoilTemperatureMod.F90

                !determine the type of pedotransfer function to be used based on soil order
                !I will use the following implementation to further explore the ET problem, now
                !I set soil order to 0 for all soils. Jinyun Tang, Mar 20, 2014

                ipedof=get_ipedof(0)
                call pedotransf(ipedof, sand, clay, &
                     soilstate_inst%watsat_col(c,lev), soilstate_inst%bsw_col(c,lev), soilstate_inst%sucsat_col(c,lev), xksat)

                om_watsat         = max(0.93_r8 - 0.1_r8   *(zsoi(lev)/zsapric), 0.83_r8)
                om_b              = min(2.7_r8  + 9.3_r8   *(zsoi(lev)/zsapric), 12.0_r8)
                om_sucsat         = min(10.3_r8 - 0.2_r8   *(zsoi(lev)/zsapric), 10.1_r8)
                om_hksat          = max(0.28_r8 - 0.2799_r8*(zsoi(lev)/zsapric), xksat)

                soilstate_inst%bd_col(c,lev)        = (1._r8 - soilstate_inst%watsat_col(c,lev))*params_inst%pd
                ! do not allow watsat_sf to push watsat above 0.93
                soilstate_inst%watsat_col(c,lev)    = min(params_inst%watsat_sf * ( (1._r8 - om_frac) * &
                                                      soilstate_inst%watsat_col(c,lev) + om_watsat*om_frac ), 0.93_r8)
                tkm                                 = (1._r8-om_frac) * (params_inst%tkd_sand*sand+params_inst%tkd_clay*clay)/ &
                                                      (sand+clay)+params_inst%tkm_om*om_frac ! W/(m K)
                soilstate_inst%bsw_col(c,lev)       = params_inst%bsw_sf * ( (1._r8-om_frac) * &
                                                      (2.91_r8 + 0.159_r8*clay) + om_frac*om_b )
                soilstate_inst%sucsat_col(c,lev)    = params_inst%sucsat_sf * ( (1._r8-om_frac) * &
                                                      soilstate_inst%sucsat_col(c,lev) + om_sucsat*om_frac ) 
                soilstate_inst%hksat_min_col(c,lev) = xksat

                ! perc_frac is zero unless perf_frac greater than percolation threshold
                if (om_frac > pcalpha) then
                   perc_norm=(1._r8 - pcalpha)**(-pcbeta)
                   perc_frac=perc_norm*(om_frac - pcalpha)**pcbeta
                else
                   perc_frac=0._r8
                endif

                ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
                uncon_frac=(1._r8-om_frac)+(1._r8-perc_frac)*om_frac

                ! uncon_hksat is series addition of mineral/organic conductivites
                if (om_frac < 1._r8) then
                   uncon_hksat=uncon_frac/((1._r8-om_frac)/xksat &
                        +((1._r8-perc_frac)*om_frac)/om_hksat)
                else
                   uncon_hksat = 0._r8
                end if
                soilstate_inst%hksat_col(c,lev)  = params_inst%hksat_sf * ( uncon_frac*uncon_hksat + &
                                                   (perc_frac*om_frac)*om_hksat )

                soilstate_inst%tkmg_col(c,lev)   = tkm ** (1._r8- soilstate_inst%watsat_col(c,lev))           

                soilstate_inst%tksatu_col(c,lev) = soilstate_inst%tkmg_col(c,lev)*0.57_r8**soilstate_inst%watsat_col(c,lev)

                soilstate_inst%tkdry_col(c,lev)  = ((0.135_r8*soilstate_inst%bd_col(c,lev) + 64.7_r8) / &
                     (params_inst%pd - 0.947_r8*soilstate_inst%bd_col(c,lev)))*(1._r8-om_frac) + params_inst%tkd_om*om_frac  

                soilstate_inst%csol_col(c,lev)   = ((1._r8-om_frac)*(params_inst%csol_sand*sand+ &
                     params_inst%csol_clay*clay) / (sand+clay) + params_inst%csol_om*om_frac)*1.e6_r8  ! J/(m3 K)

                soilstate_inst%watdry_col(c,lev) = soilstate_inst%watsat_col(c,lev) * &
                     (316230._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 
                soilstate_inst%watopt_col(c,lev) = soilstate_inst%watsat_col(c,lev) * &
                     (158490._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

                !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
                ! water content at field capacity, defined as hk = 0.1 mm/day
                ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / secspday (day/sec)
                soilstate_inst%watfc_col(c,lev) = soilstate_inst%watsat_col(c,lev) * &
                     (0.1_r8 / (soilstate_inst%hksat_col(c,lev)*secspday))**(1._r8/(2._r8*soilstate_inst%bsw_col(c,lev)+3._r8))
             end if
          end do

          ! Urban pervious and impervious road
          if (col%itype(c) == icol_road_imperv) then
             ! Impervious road layers -- same as above except set watdry and watopt as missing
             do lev = 1,nlevgrnd
                soilstate_inst%watdry_col(c,lev) = spval 
                soilstate_inst%watopt_col(c,lev) = spval 
             end do
          else if (col%itype(c) == icol_road_perv) then 
             ! pervious road layers  - set in UrbanInitTimeConst
          end if

       end if
    end do

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: lake
    ! --------------------------------------------------------------------

    do c = begc, endc
       g = col%gridcell(c)
       l = col%landunit(c)

       if (lun%itype(l)==istdlak) then

          do lev = 1,nlevgrnd
             if ( lev <= nlevsoi )then
                clay    =  soilstate_inst%cellclay_col(c,lev)
                sand    =  soilstate_inst%cellsand_col(c,lev)
                if ( organic_frac_squared )then
                   om_frac = min( params_inst%om_frac_sf*(soilstate_inst%cellorg_col(c,lev)/organic_max)**2._r8, 1._r8)
                else
                   om_frac = min(params_inst%om_frac_sf*soilstate_inst%cellorg_col(c,lev)/organic_max, 1._r8)
                end if
             else
                clay    = soilstate_inst%cellclay_col(c,nlevsoi)
                sand    = soilstate_inst%cellsand_col(c,nlevsoi)
                om_frac = 0.0_r8
             end if

             soilstate_inst%watsat_col(c,lev) = 0.489_r8 - 0.00126_r8*sand

             soilstate_inst%bsw_col(c,lev)    = 2.91 + 0.159*clay

             soilstate_inst%sucsat_col(c,lev) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )

             bd = (1._r8-soilstate_inst%watsat_col(c,lev))*params_inst%pd

             ! do not allow watsat_sf to push watsat above 0.93
             soilstate_inst%watsat_col(c,lev) = min(params_inst%watsat_sf * ( (1._r8 - om_frac) * &
                    soilstate_inst%watsat_col(c,lev) + om_watsat_lake * om_frac), 0.93_r8)

             tkm = (1._r8-om_frac)*(params_inst%tkd_sand*sand+params_inst%tkd_clay*clay)/(sand+clay) + &
                   params_inst%tkm_om * om_frac ! W/(m K)

             soilstate_inst%bsw_col(c,lev)    = params_inst%bsw_sf * ( (1._r8-om_frac) * &
                   (2.91_r8 + 0.159_r8*clay) + om_frac * om_b_lake )

             soilstate_inst%sucsat_col(c,lev) = params_inst%sucsat_sf * ( (1._r8-om_frac) * &
                   soilstate_inst%sucsat_col(c,lev) + om_sucsat_lake * om_frac )

             xksat = 0.0070556_r8 *( 10._r8**(-0.884_r8+0.0153_r8*sand) ) ! mm/s

             ! perc_frac is zero unless perf_frac greater than percolation threshold
             if (om_frac > pc_lake) then
                perc_norm = (1._r8 - pc_lake)**(-pcbeta)
                perc_frac = perc_norm*(om_frac - pc_lake)**pcbeta
             else
                perc_frac = 0._r8
             endif

             ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
             uncon_frac = (1._r8-om_frac) + (1._r8-perc_frac)*om_frac

             ! uncon_hksat is series addition of mineral/organic conductivites
             if (om_frac < 1._r8) then
                xksat = 0.0070556_r8 *( 10._r8**(-0.884_r8+0.0153_r8*sand) ) ! mm/s
                uncon_hksat = uncon_frac/((1._r8-om_frac)/xksat + ((1._r8-perc_frac)*om_frac)/om_hksat_lake)
             else
                uncon_hksat = 0._r8
             end if

             soilstate_inst%hksat_col(c,lev)  = params_inst%hksat_sf * ( uncon_frac*uncon_hksat + &
                                       (perc_frac*om_frac)*om_hksat_lake )
             soilstate_inst%tkmg_col(c,lev)   = tkm ** (1._r8- soilstate_inst%watsat_col(c,lev))
             soilstate_inst%tksatu_col(c,lev) = soilstate_inst%tkmg_col(c,lev)*0.57_r8**soilstate_inst%watsat_col(c,lev)
             soilstate_inst%tkdry_col(c,lev)  = ((0.135_r8*bd + 64.7_r8) / (params_inst%pd - 0.947_r8*bd))*(1._r8-om_frac) + &
                                       params_inst%tkd_om * om_frac
             soilstate_inst%csol_col(c,lev)   = ((1._r8-om_frac)*(params_inst%csol_sand*sand+ &
                                       params_inst%csol_clay*clay) / (sand+clay) + params_inst%csol_om * om_frac)*1.e6_r8 ! J/(m3 K)
             soilstate_inst%watdry_col(c,lev) = soilstate_inst%watsat_col(c,lev) &
                  * (316230._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev))
             soilstate_inst%watopt_col(c,lev) = soilstate_inst%watsat_col(c,lev) &
                  * (158490._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev))

             !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
             ! water content at field capacity, defined as hk = 0.1 mm/day
             ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / (# seconds/day)
             soilstate_inst%watfc_col(c,lev) = soilstate_inst%watsat_col(c,lev) * (0.1_r8 / &
                  (soilstate_inst%hksat_col(c,lev)*secspday))**(1._r8/(2._r8*soilstate_inst%bsw_col(c,lev)+3._r8))
          end do
       endif

    end do

    ! --------------------------------------------------------------------
    ! Initialize threshold soil moisture, and mass fraction of clay as
    ! scaling coefficient of dust emission flux (kg/m2/s) in each DustEmisType
    ! module. See the comments in each function.
    ! Zender suggested that threshold soil moisture is tunable (see comment
    ! inside ThresholdSoilMoistZender2003). dmleung further add dust_moist_fact
    ! ofr modelers to tune the threshold soil moisture. The resulting tuning 
    ! factor is thus a = dust_moist_fact / (clay3d). dmleung 30 Sep 2024
    ! --------------------------------------------------------------------

    do c = begc,endc
       g = col%gridcell(c)

       !soilstate_inst%gwc_thr_col(c) = ThresholdSoilMoistZender2003( clay3d(g,1) )
       if ( is_dust_emis_leung() )then
          soilstate_inst%mss_frc_cly_vld_col(c) = MassFracClayLeung2023( clay3d(g,1) )
          dust_moist_fact = 1.0_r8   ! change this into a namelist variable later., currrently not used but could be in the future
       else
          soilstate_inst%mss_frc_cly_vld_col(c) = MassFracClay( clay3d(g,1) )
          dust_moist_fact = 1.0_r8
       end if
       soilstate_inst%gwc_thr_col(c) = dust_moist_fact * ThresholdSoilMoistZender2003( clay3d(g,1) )

    end do

    ! --------------------------------------------------------------------
    ! Deallocate memory
    ! --------------------------------------------------------------------

    deallocate(sand3d, clay3d, organic3d)
    deallocate(zisoifl, zsoifl)

  end subroutine SoilStateInitTimeConst

  !------------------------------------------------------------------------------

  real(r8) function ThresholdSoilMoistZender2003( clay )
  !------------------------------------------------------------------------------
  !
  ! Calculate the threshold gravimetric water content needed for dust emission, based on clay content
  ! This was the original equation with a = 1 / (%clay) being the tuning factor for soil
  ! moisture effect in Zender's 2003 dust emission scheme (only for top layer).
  ! dmleung further added dust_moist_fact for more flexibility in tuning, so the tuning factor here
  ! is a = dust_moist_fact / (%clay). dmleung added dust_moist_fact on 30 Sep 2024.
  !
  ! 0.17 and 0.14 are fitting coefficients in Fecan et al. (1999), and 0.01 is used to
  ! convert surface clay from percentage to fraction.
  ! The equation comes from Eq. 14 of Fecan et al. (1999; https://doi.org/10.1007/s00585-999-0149-7).
  !
  ! NOTE: dmleung 19 Feb 2024.
  !------------------------------------------------------------------------------
  ! For future developments Danny M. Leung decided (Dec, 2023) that the Leung et al. (2023) o
  ! dust emission scheme in the CESM will use Zender's tuning as well, which overall
  ! encourages more dust emissions from seriamid and more marginal dust sources.
  ! Another advantage of using this tuning factor instead of a = 1 is that the dust emission
  ! threshold is linearly dependent on the clay fraction instead of parabolically dependent
  ! on clay fraction as in the above line. This means that dust emission becomes a little
  ! less sensitive to clay content (soil texture).
  !
  ! Also see the notes below for ThresholdSoilMoistKok2014.
  !
  ! NOTE on Leung dust emissions: dmleung Jul 5 2024:
  !
  ! dmleung followed Zender (2003) DUST scheme to again avoid dust flux being too sensitive to the choice
  ! of clay dataset. This is different from what Leung et al. (2023) intended to do.
  ! NOTE: This might need to be adjusted for tuning in the future.
  !
  !------------------------------------------------------------------------------
      real(r8), intent(IN) :: clay ! Fraction of clay in the soil (%)

      if ( clay < 0.0_r8 .or. clay > 100.0_r8 )then
         ThresholdSoilMoistZender2003 = nan
         call endrun( 'Clay fraction is out of bounds (0 to 100)')
         return
      end if
      ThresholdSoilMoistZender2003 = 0.17_r8 + 0.14_r8 * clay * 0.01_r8
  end function ThresholdSoilMoistZender2003

  !------------------------------------------------------------------------------

  real(r8) function ThresholdSoilMoistKok2014( clay )
  !------------------------------------------------------------------------------
  ! Calculate the threshold soil moisture needed for dust emission, based on clay content
  !
  ! NOTE: dmleung 24 May 2024.
  !
  ! The below calculates the threshold gravimetric water content for the dust emission
  ! calculation in DustEmis. The equation comes from Eq. 14 of Fecan et al.
  ! (1999; https://doi.org/10.1007/s00585-999-0149-7).
  ! gwc_thr_col = 0.17*clay3d + 0.0014*(clay3d**2), and we only concern the topmost
  ! soil layer.  Charlie Zender later on added a tuning factor (a) such that the
  ! equation becomes gwc_thr_col = a*[0.17*clay3d + 0.0014*(clay3d**2)].
  ! (Zender et al., 2003a; https://doi.org/10.1029/2002JD002775)
  ! Kok et al. (2014a, b) chose to use a = 1. Resulting in this function
  ! Charlie Zender (2003a) chose:  a = 1/clay3d, which gives the ThresholdSoilMoistZender2003
  ! function above.
  !
  !------------------------------------------------------------------------------
      real(r8), intent(IN) :: clay ! Fraction of clay in the soil (%)

      ThresholdSoilMoistKok2014 = 0.01_r8*(0.17_r8*clay + 0.0014_r8*clay*clay)
  end function ThresholdSoilMoistKok2014

  !------------------------------------------------------------------------------

  real(r8) function MassFracClay( clay )
  ! Calculate the mass fraction of clay needed for dust emission, based on clay content
      real(r8), intent(IN) :: clay ! Fraction of lay in the soil (%)

      MassFracClay = min(clay * 0.01_r8, 0.20_r8)
  end function MassFracClay

  !------------------------------------------------------------------------------

  real(r8) function MassFracClayLeung2023( clay )
  ! Calculate the mass fraction of clay needed for dust emission, based on clay content
  ! Based on the base Zender_2003 version, with a slight modification for Leung_2023
  ! dmleung modified 5 Jul 2024, reducing sensitivity of dust emission
  ! flux to clay fraction.
  ! NOTE: This might need to be adjusted for tuning in the future.
      real(r8), intent(IN) :: clay ! Fraction of lay in the soil (%)

      MassFracClayLeung2023 = 0.1_r8 + MassFracClay( clay ) * 0.1_r8 / 0.20_r8   ! dmleung added this line to reduce the sensitivity of dust emission flux to clay fraction in DUSTMod. 5 Jul 2024
  end function MassFracClayLeung2023

  !------------------------------------------------------------------------------

end module SoilStateInitTimeConstMod
