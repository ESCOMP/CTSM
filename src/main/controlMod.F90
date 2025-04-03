module controlMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module which initializes run control variables. The following possible
  ! namelist variables are set default values and possibly read in on startup
  !
  ! Note: For definitions of namelist variables see
  !       ../../bld/namelist_files/namelist_definition.xml
  !       Display the file in a browser to see it neatly formatted in html.
  !
  ! !USES:
  use shr_kind_mod                     , only: r8 => shr_kind_r8
  use shr_nl_mod                       , only: shr_nl_find_group_name
  use shr_const_mod                    , only: SHR_CONST_CDAY
  use shr_log_mod                      , only: errMsg => shr_log_errMsg
  use abortutils                       , only: endrun
  use spmdMod                          , only: masterproc, mpicom
  use spmdMod                          , only: MPI_CHARACTER, MPI_INTEGER, MPI_LOGICAL, MPI_REAL8
  use decompInitMod                    , only: clump_pproc
  use clm_varcon                       , only: h2osno_max
  use clm_varpar                       , only: maxpatch_glc, numrad, nlevsno
  use fileutils                        , only: getavu, relavu, get_filename
  use histFileMod                      , only: max_tapes, max_namlen
  use histFileMod                      , only: hist_empty_htapes, hist_dov2xy, hist_avgflag_pertape, hist_type1d_pertape
  use histFileMod                      , only: hist_nhtfrq, hist_ndens, hist_mfilt, hist_fincl1, hist_fincl2, hist_fincl3
  use histFileMod                      , only: hist_fincl4, hist_fincl5, hist_fincl6, hist_fincl7, hist_fincl8
  use histFileMod                      , only: hist_fincl9, hist_fincl10
  use histFileMod                      , only: hist_fexcl1, hist_fexcl2, hist_fexcl3,  hist_fexcl4, hist_fexcl5, hist_fexcl6
  use histFileMod                      , only: hist_fexcl7, hist_fexcl8, hist_fexcl9, hist_fexcl10
  use initInterpMod                    , only: initInterp_readnl
  use LakeCon                          , only: deepmixing_depthcrit, deepmixing_mixfact
  use CanopyfluxesMod                  , only: perchroot, perchroot_alt
  use CanopyHydrologyMod               , only: CanopyHydrology_readnl
  use SurfaceAlbedoMod                 , only: SurfaceAlbedo_readnl
  use SurfaceResistanceMod             , only: soil_resistance_readNL
  use SnowHydrologyMod                 , only: SnowHydrology_readnl
  use SurfaceAlbedoMod                 , only: albice, lake_melt_icealb
  use UrbanParamsType                  , only: UrbanReadNML
  use HumanIndexMod                    , only: HumanIndexReadNML
  use CNPrecisionControlMod            , only: CNPrecisionControlReadNML
  use CNSharedParamsMod                , only: use_fun, use_matrixcn
  use CIsoAtmTimeseriesMod             , only: use_c14_bombspike, atm_c14_filename, use_c13_timeseries, atm_c13_filename
  use SoilBiogeochemDecompCascadeConType, only : use_soil_matrixcn
  use SoilBiogeochemCompetitionMod     , only: suplnitro, suplnNon
  use SoilBiogeochemLittVertTranspMod  , only: som_adv_flux, max_depth_cryoturb
  use SoilBiogeochemVerticalProfileMod , only: surfprof_exp
  use SoilBiogeochemNitrifDenitrifMod  , only: no_frozen_nitrif_denitrif
  use SoilHydrologyMod                 , only: soilHydReadNML, hillslope_hydrology_ReadNML
  use CNFireFactoryMod                 , only: CNFireReadNML
  use CanopyFluxesMod                  , only: CanopyFluxesReadNML
  use shr_drydep_mod                   , only: n_drydep
  use clm_varctl
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: control_setNL ! Set namelist filename
  public :: control_init  ! initial run control information
  public :: control_print ! print run control information
  public :: check_missing_initdata_status  ! check for missing finidat_interp_dest .status file
  !
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: apply_use_init_interp  ! apply the use_init_interp namelist option, if set
  !
  ! !PRIVATE TYPES:
  character(len=  7) :: runtyp(4)                        ! run type
  character(len=fname_len) :: NLFilename = 'lnd.stdin' ! Namelist filename

#if (defined _OPENMP)
  integer, external :: omp_get_max_threads  ! max number of threads that can execute concurrently in a single parallel region
#endif

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine control_setNL( NLfile )
    !
    ! !DESCRIPTION:
    ! Set the namelist filename to use
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFile ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'control_setNL'  ! subroutine name
    logical :: lexist                               ! File exists
    !------------------------------------------------------------------------

    ! Error checking...
    if ( len_trim(NLFile) == 0 )then
       call endrun(msg=' error: nlfilename entered is not set'//errMsg(sourcefile, __LINE__))
    end if
    inquire (file = trim(NLFile), exist = lexist)
    if ( .not. lexist )then
       call endrun(msg=' error: NLfilename entered does NOT exist:'//&
            trim(NLFile)//errMsg(sourcefile, __LINE__))
    end if
    if ( len_trim(NLFile) > len(NLFilename) )then
       call endrun(msg=' error: entered NLFile is too long'//errMsg(sourcefile, __LINE__))
    end if
    ! Set the filename
    NLFilename = NLFile
    NLFilename_in = NLFilename   ! For use in external namelists and to avoid creating dependencies on controlMod
  end subroutine control_setNL

  !------------------------------------------------------------------------

  subroutine control_init(dtime)
    !
    ! !DESCRIPTION:
    ! Initialize CLM run control information
    !
    ! !USES:
    use CNMRespMod                       , only : CNMRespReadNML
    use CNFUNMod                         , only : CNFUNReadNML
    use CNNDynamicsMod                   , only : CNNDynamicsReadNML
    use CNPhenologyMod                   , only : CNPhenologyReadNML
    use landunit_varcon                  , only : max_lunit
    use CNSoilMatrixMod                  , only : CNSoilMatrixInit
    !
    ! ARGUMENTS
    integer, intent(in) :: dtime    ! model time step (seconds)

    ! !LOCAL VARIABLES:
    integer :: i                    ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    logical :: use_init_interp      ! Apply initInterp to the file given by finidat
    !------------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! CLM namelist settings

    namelist /clm_inparm/ &
         fatmlndfrc, finidat, nrevsn, &
         use_init_interp, compname

    ! Input datasets

    namelist /clm_inparm/  &
         fsurdat, hillslope_file, &
         paramfile, fsnowoptics, fsnowaging

    ! History, restart options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fincl7,  hist_fincl8,              &
         hist_fincl9,  hist_fincl10,             &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6, &
         hist_fexcl7,  hist_fexcl8,              &
         hist_fexcl9,  hist_fexcl10
    namelist /clm_inparm/ hist_wrtch4diag, hist_fields_list_file

    ! BGC info

    namelist /clm_inparm/  &
         suplnitro
    namelist /clm_inparm/ &
         nfix_timeconst
    namelist /clm_inparm/ &
         spinup_state, override_bgc_restart_mismatch_dump

    namelist /clm_inparm / &
         co2_type

    namelist /clm_inparm / &
         perchroot, perchroot_alt

    namelist /clm_inparm / &
         anoxia, use_fun

    namelist /clm_inparm / &
         deepmixing_depthcrit, deepmixing_mixfact, lake_melt_icealb

    ! CN Matrix solution
    namelist /clm_inparm / &
         use_matrixcn, use_soil_matrixcn, hist_wrt_matrixcn_diag, spinup_matrixcn, nyr_forcing, nyr_sasu, iloop_avg

    ! lake_melt_icealb is of dimension numrad

    ! Glacier_mec info
    namelist /clm_inparm/ &
         maxpatch_glc, glc_do_dynglacier, &
         glc_snow_persistence_max_days, &
         nlevsno, h2osno_max

    ! Other options

    namelist /clm_inparm/  &
         clump_pproc, &
         create_crop_landunit, nsegspc, co2_ppmv, &
         albice, soil_layerstruct_predefined, soil_layerstruct_userdefined, &
         soil_layerstruct_userdefined_nlevsoi, use_subgrid_fluxes, &
         snow_thermal_cond_method, snow_cover_fraction_method, &
         irrigate, run_zero_weight_urban, all_active, &
         crop_fsat_equals_zero, for_testing_run_ncdiopio_tests, &
         for_testing_use_second_grain_pool, for_testing_use_repr_structure_pool, &
         for_testing_no_crop_seed_replenishment, &
         z0param_method, use_z0m_snowmelt

    ! vertical soil mixing variables
    namelist /clm_inparm/  &
         som_adv_flux, max_depth_cryoturb

    ! C and N input vertical profiles
    namelist /clm_inparm/  &
          surfprof_exp

    namelist /clm_inparm/ no_frozen_nitrif_denitrif

    namelist /clm_inparm/ use_c13, use_c14, for_testing_allow_interp_non_ciso_to_ciso


    ! FATES Flags
    namelist /clm_inparm/ fates_paramfile, use_fates,   &
          fates_spitfire_mode, fates_harvest_mode,      &
          use_fates_planthydro, use_fates_ed_st3,       &
          use_fates_cohort_age_tracking,                &
          use_fates_ed_prescribed_phys,                 &
          use_fates_inventory_init,                     &
          use_fates_fixed_biogeog,                      &
          use_fates_nocomp,                             &
          use_fates_sp,                                 &
          use_fates_luh,                                &
          use_fates_lupft,                              &
          use_fates_potentialveg,                       &
          fluh_timeseries,                              &
          flandusepftdat,                               &
          fates_inventory_ctrl_filename,                &
          fates_stomatal_model,                         &
          fates_stomatal_assimilation,                  &
          fates_leafresp_model,                         &
          fates_cstarvation_model,                      &
          fates_regeneration_model,                     &
          fates_radiation_model,                        &
          fates_hydro_solver,                           &
          fates_electron_transport_model,               &
          fates_parteh_mode,                            &
          fates_seeddisp_cadence,                       &
          use_fates_tree_damage,                        &
          use_fates_daylength_factor,                   &
          fates_photosynth_acclimation,                 &
          fates_history_dimlevel

    ! Ozone vegetation stress method
    namelist / clm_inparm / o3_veg_stress_method

    ! CLM 5.0 nitrogen flags
    namelist /clm_inparm/ use_flexibleCN, use_luna

    namelist /clm_nitrogen/ MM_Nuptake_opt, &
         CNratio_floating, lnc_opt, reduce_dayl_factor, vcmax_opt, &
         CN_evergreen_phenology_opt, carbon_resp_opt

    namelist /clm_inparm/ use_soil_moisture_streams

    ! excess ice flag
    namelist /clm_inparm/ use_excess_ice 

    namelist /clm_inparm/ use_lai_streams

    namelist /clm_inparm/ use_bedrock

    namelist /clm_inparm/ use_biomass_heat_storage

    namelist /clm_inparm/ use_hillslope

    namelist /clm_inparm/ downscale_hillslope_meteorology

    namelist /clm_inparm/ use_hillslope_routing

    namelist /clm_inparm/ hillslope_fsat_equals_zero

    namelist /clm_inparm/ use_hydrstress

    namelist /clm_inparm/  &
         use_c14_bombspike, atm_c14_filename, use_c13_timeseries, atm_c13_filename

    ! All old cpp-ifdefs are below and have been converted to namelist variables

    namelist /clm_inparm/ convert_ocean_to_land

    ! Number of dominant pfts and landunits. Enhance ctsm performance by
    ! reducing the number of active pfts to n_dom_pfts and
    ! active landunits to n_dom_landunits.
    ! Also choose to collapse the urban landunits to the dominant urban
    ! landunit by setting collapse_urban = .true.
    namelist /clm_inparm/ n_dom_pfts
    namelist /clm_inparm/ n_dom_landunits
    namelist /clm_inparm/ collapse_urban

    ! Thresholds above which the model keeps the soil, crop, glacier, lake,
    ! wetland, and urban landunits
    namelist /clm_inparm/ toosmall_soil, toosmall_crop, toosmall_glacier
    namelist /clm_inparm/ toosmall_lake, toosmall_wetland, toosmall_urban

    ! flag for SSRE diagnostic
    namelist /clm_inparm/ use_SSRE

    namelist /clm_inparm/ &
         use_lch4, use_nitrif_denitrif, use_extralakelayers, &
         use_vichydro, use_cn, use_cndv, use_crop, use_fertilizer, &
         use_grainproduct, use_snicar_frc, use_vancouver, use_mexicocity, use_noio, &
         use_nguardrail, crop_residue_removal_frac, flush_gdd20

    ! SNICAR
    namelist /clm_inparm/ &
         snicar_numrad_snw, snicar_solarspec, snicar_dust_optics, &
         snicar_use_aerosol, snicar_snw_shape, snicar_snobc_intmix, &
         snicar_snodst_intmix, do_sno_oc

    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize run control settings .....'
    endif

    finidat_interp_dest = './init_generated_files/finidat_interp_dest'//trim(inst_suffix)//'.nc'
    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'
    runtyp(nsrBranch   + 1) = 'branch '

    if(use_fates)then
       use_matrixcn = .false.
       use_soil_matrixcn = .false.
       hist_wrt_matrixcn_diag = .false.
       spinup_matrixcn = .false.
    end if
    nyr_forcing = 10

    ! Set clumps per procoessor

#if (defined _OPENMP)
    clump_pproc = omp_get_max_threads()
#else
    clump_pproc = 1
#endif

    use_init_interp = .false.
    if (masterproc) then

       ! ----------------------------------------------------------------------
       ! Read namelist from standard input.
       ! ----------------------------------------------------------------------

       if ( len_trim(NLFilename) == 0  )then
          call endrun(msg=' error: nlfilename not set'//errMsg(sourcefile, __LINE__))
       end if
       unitn = getavu()
       write(iulog,*) 'Read in clm_inparm namelist from: ', trim(NLFilename)
       open( unitn, file=trim(NLFilename), status='old' )
       call shr_nl_find_group_name(unitn, 'clm_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='ERROR reading clm_inparm namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='ERROR finding clm_inparm namelist'//errMsg(sourcefile, __LINE__))
       end if
       call shr_nl_find_group_name(unitn, 'clm_nitrogen', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_nitrogen, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='ERROR reading clm_nitrogen namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='ERROR finding clm_nitrogen namelist'//errMsg(sourcefile, __LINE__))
       end if

       call relavu( unitn )

       ! ----------------------------------------------------------------------
       ! Process some namelist variables, and perform consistency checks
       ! ----------------------------------------------------------------------

       ! History and restart files (dependent on settings of dtime)
       do i = 1, max_tapes
          if (hist_nhtfrq(i) < 0) then
             hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*dtime))
          endif
       end do

       if (use_init_interp) then
          call apply_use_init_interp(finidat_interp_dest, finidat, finidat_interp_source)
       end if

       if (maxpatch_glc <= 0) then
          call endrun(msg=' ERROR: maxpatch_glc must be at least 1 ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (n_dom_pfts < 0) then
          call endrun(msg=' ERROR: expecting n_dom_pfts between 0 and 14 where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if
       if (n_dom_landunits < 0 .or. n_dom_landunits > max_lunit) then
          call endrun(msg=' ERROR: expecting n_dom_landunits between 0 and  max_lunit where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (toosmall_soil < 0._r8 .or. toosmall_soil > 100._r8) then
          call endrun(msg=' ERROR: expecting toosmall_soil between 0._r8 and 100._r8 where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (toosmall_crop < 0._r8 .or. toosmall_crop > 100._r8) then
          call endrun(msg=' ERROR: expecting toosmall_crop between 0._r8 and 100._r8 where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (toosmall_glacier < 0._r8 .or. toosmall_glacier > 100._r8) then
          call endrun(msg=' ERROR: expecting toosmall_glacier between 0._r8 and 100._r8 where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (toosmall_lake < 0._r8 .or. toosmall_lake > 100._r8) then
          call endrun(msg=' ERROR: expecting toosmall_lake between 0._r8 and 100._r8 where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (toosmall_wetland < 0._r8 .or. toosmall_wetland > 100._r8) then
          call endrun(msg=' ERROR: expecting toosmall_wetland between 0._r8 and 100._r8 where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (toosmall_urban < 0._r8 .or. toosmall_urban > 100._r8) then
          call endrun(msg=' ERROR: expecting toosmall_urban between 0._r8 and 100._r8 where 0 is the default value that tells the model to do nothing ' // &
               errMsg(sourcefile, __LINE__))
       end if

       if (glc_do_dynglacier) then
          if (collapse_urban) then
             call endrun(msg='ERROR: glc_do_dynglacier is incompatible &
                              with collapse_urban = .true.' // &
                              errMsg(sourcefile, __LINE__))
          end if
          if (n_dom_pfts > 0 .or. n_dom_landunits > 0 &
              .or. toosmall_soil > 0._r8 .or. toosmall_crop > 0._r8 &
              .or. toosmall_glacier > 0._r8 .or. toosmall_lake > 0._r8 &
              .or. toosmall_wetland > 0._r8 .or. toosmall_urban > 0._r8) then
             call endrun(msg='ERROR: glc_do_dynglacier is incompatible &
                              with any of the following set to > 0: &
                              n_dom_pfts > 0, n_dom_landunits > 0, &
                              toosmall_soil > 0._r8, toosmall_crop > 0._r8, &
                              toosmall_glacier > 0._r8, toosmall_lake > 0._r8, &
                              toosmall_wetland > 0._r8, toosmall_urban > 0._r8.' // &
                              errMsg(sourcefile, __LINE__))
          end if
       end if

       if (use_crop .and. .not. create_crop_landunit) then
          call endrun(msg=' ERROR: prognostic crop Patches require create_crop_landunit=.true.'//&
            errMsg(sourcefile, __LINE__))
       end if

       if (use_lch4 ) then
          anoxia = .true.
       else
          anoxia = .false.
       end if

       ! ----------------------------------------------------------------------
       ! Check compatibility with the FATES model
       if ( use_fates ) then

          if(use_fates_sp) then
             use_fates_bgc = .false.
          else
             use_fates_bgc = .true.
          end if
          
          if (fates_parteh_mode == 1 .and. suplnitro == suplnNon .and. use_fates_bgc )then
             write(iulog,*) ' When FATES with fates_parteh_mode == 1 (ie carbon only mode),'
             write(iulog,*) '  you must have supplemental nitrogen turned on, there will be'
             write(iulog,*) '  no nitrogen dynamics with the plants, and therefore no'
             write(iulog,*) '  meaningful limitations to nitrogen.'
             call endrun(msg=' ERROR: fates_parteh_mode=1 must have suplnitro set to suplnAll.'//&
                   errMsg(sourcefile, __LINE__))
          end if
          
          if ( use_cn) then
             call endrun(msg=' ERROR: use_cn and use_fates cannot both be set to true.'//&
                   errMsg(sourcefile, __LINE__))
          end if

          if ( use_hydrstress) then
             call endrun(msg=' ERROR: use_hydrstress and use_fates cannot both be set to true.'//&
                   errMsg(sourcefile, __LINE__))
          end if

          if ( use_crop ) then
             call endrun(msg=' ERROR: use_crop and use_fates cannot both be set to true.'//&
                   errMsg(sourcefile, __LINE__))
          end if

          if( use_luna ) then
             call endrun(msg=' ERROR: luna is not compatible with FATES.'//&
                  errMsg(sourcefile, __LINE__))
          end if

          if (o3_veg_stress_method /= 'unset' ) then
             call endrun(msg=' ERROR: ozone is not compatible with FATES.'//&
                  errMsg(sourcefile, __LINE__))
          end if

          if (use_c13 .or. use_c14) then
             call endrun(msg=' ERROR: C13 and C14 dynamics are not compatible with FATES.'//&
                  errMsg(sourcefile, __LINE__))
          end if

          if (z0param_method == 'Meier2022') then
             call endrun(msg=' ERROR: Surface roughness parameterization Meier2022 is not compatible with FATES.'//&
                  errMsg(sourcefile, __LINE__))
          end if

       else
          
          ! These do default to false anyway, but this emphasizes they
          ! are false when fates is false
          use_fates_sp  = .false.
          use_fates_bgc = .false.
          
       end if

       ! Check compatibility with use_lai_streams
       if (use_lai_streams) then
        if ((use_fates .and. .not. use_fates_sp) .or. use_cn) then 
          call endrun(msg=' ERROR: cannot use LAI streams unless in SP mode (use_cn = .false. or use_fates_sp=.true.).'//&
                  errMsg(sourcefile, __LINE__))
        end if 
       end if 

       ! If nfix_timeconst is equal to the junk default value, then it was not specified
       ! by the user namelist and we need to assign it the correct default value. If the
       ! user specified it in the namelist, we leave it alone.

       if (nfix_timeconst == -1.2345_r8) then
          if (use_nitrif_denitrif) then
             nfix_timeconst = 10._r8
          else
             nfix_timeconst = 0._r8
          end if
       end if

       ! If nlevsno or h2osno_max are equal to their junk
       ! default value, then they were not specified by the user namelist and we generate
       ! an error message. Also check nlevsno for bounds.
       if (nlevsno < 3 .or. nlevsno > 12)  then
          write(iulog,*)'ERROR: nlevsno = ',nlevsno,' is not supported, must be in range 3-12.'
          call endrun(msg=' ERROR: invalid value for nlevsno in CLM namelist. '//&
               errMsg(sourcefile, __LINE__))
       endif
       if (h2osno_max <= 0.0_r8) then
          write(iulog,*)'ERROR: h2osno_max = ',h2osno_max,' is not supported, must be greater than 0.0.'
          call endrun(msg=' ERROR: invalid value for h2osno_max in CLM namelist. '//&
               errMsg(sourcefile, __LINE__))
       endif

    endif   ! end of if-masterproc if-block

    ! ----------------------------------------------------------------------
    ! Read in other namelists for other modules
    ! ----------------------------------------------------------------------

    call mpi_bcast (use_init_interp, 1, MPI_LOGICAL, 0, mpicom, ierr)
    if (use_init_interp) then
       call initInterp_readnl( NLFilename )
    end if

    !I call init_hydrology to set up default hydrology sub-module methods.
    !For future version, I suggest to  put the following two calls inside their
    !own modules, which are called from their own initializing methods
    call init_hydrology( NLFilename )

    call soil_resistance_readnl ( NLFilename )
    call CanopyFluxesReadNML    ( NLFilename )
    call CanopyHydrology_readnl ( NLFilename )
    call SurfaceAlbedo_readnl   ( NLFilename )
    call SnowHydrology_readnl   ( NLFilename )
    call UrbanReadNML           ( NLFilename )
    call HumanIndexReadNML      ( NLFilename )

    ! ----------------------------------------------------------------------
    ! Broadcast all control information if appropriate
    ! ----------------------------------------------------------------------

    call control_spmd()

    ! ----------------------------------------------------------------------
    ! Read in other namelists that are dependent on other namelist setttings
    ! ----------------------------------------------------------------------

    if ( use_fun ) then
       call CNMRespReadNML( NLFilename )
       call CNFUNReadNML( NLFilename )
    end if

    call soilHydReadNML(   NLFilename )
    if ( use_hillslope ) then
       call hillslope_hydrology_ReadNML(   NLFilename )
    endif
    if ( use_cn ) then
       call CNFireReadNML(             NLFilename )
       call CNPrecisionControlReadNML( NLFilename )
       call CNNDynamicsReadNML       ( NLFilename )
       call CNPhenologyReadNML       ( NLFilename )
    end if

    ! CN soil matrix

    call CNSoilMatrixInit()

    ! ----------------------------------------------------------------------
    ! consistency checks
    ! ----------------------------------------------------------------------

    ! Consistency settings for co2 type
    if (co2_type /= 'constant' .and. co2_type /= 'prognostic' .and. co2_type /= 'diagnostic') then
       write(iulog,*)'co2_type = ',co2_type,' is not supported'
       call endrun(msg=' ERROR:: choices are constant, prognostic or diagnostic'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! Check on run type
    if (nsrest == iundef) then
       call endrun(msg=' ERROR:: must set nsrest'//&
            errMsg(sourcefile, __LINE__))
    end if
    if (nsrest == nsrBranch .and. nrevsn == ' ') then
       call endrun(msg=' ERROR: need to set restart data file name'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! Consistency settings for co2_ppvm
    if ( (co2_ppmv <= 0.0_r8) .or. (co2_ppmv > 3000.0_r8) ) then
       call endrun(msg=' ERROR: co2_ppmv is out of a reasonable range'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! check on SNICAR BC-snow and dust-snow internal mixing
    if ( snicar_snobc_intmix .and. snicar_snodst_intmix ) then
       call endrun(msg=' ERROR: currently dust-snow and BC-snow internal mixing cannot be activated together'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! other SNICAR warnings
    if ((snicar_snw_shape /= 'sphere' .and. snicar_snw_shape /= 'hexagonal_plate') .or.  &
         snicar_solarspec /= 'mid_latitude_winter' .or.  &
         snicar_dust_optics /= 'sahara' .or.  &
         snicar_numrad_snw /= 5 .or.  &
         snicar_snodst_intmix .or.  &
         .not. snicar_use_aerosol .or.  &
         do_sno_oc) then
       call endrun(msg=' ERROR: You have selected an option that is EXPERIMENTAL, UNSUPPORTED, and UNTESTED. For guidance see namelist_defaults_ctsm.xml'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! Consistency settings for nrevsn

    if (nsrest == nsrStartup ) nrevsn = ' '
    if (nsrest == nsrContinue) nrevsn = 'set by restart pointer file file'
    if (nsrest /= nsrStartup .and. nsrest /= nsrContinue .and. nsrest /= nsrBranch ) then
       call endrun(msg=' ERROR: nsrest NOT set to a valid value'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! Single Column
    if ( single_column .and. (scmlat == rundef  .or. scmlon == rundef ) ) then
       call endrun(msg=' ERROR:: single column mode on -- but scmlat and scmlon are NOT set'//&
            errMsg(sourcefile, __LINE__))
       if (.not. use_lch4 .and. anoxia) then
          call endrun(msg='ERROR:: anoxia is turned on, but this currently requires turning on the CH4 submodel'//&
            errMsg(sourcefile, __LINE__))
       end if
    end if

    if (masterproc) then
       write(iulog,*) 'Successfully initialized run control settings'
       write(iulog,*)
    endif

  end subroutine control_init

  !------------------------------------------------------------------------
  subroutine control_spmd()
    !
    ! !DESCRIPTION:
    ! Distribute namelist data all processors. All program i/o is
    ! funnelled through the master processor. Processor 0 either
    ! reads restart/history data from the disk and distributes
    ! it to all processors, or collects data from
    ! all processors and writes it to disk.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer ier       !error code
    !-----------------------------------------------------------------------

    ! run control variables
    call mpi_bcast (caseid, len(caseid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle, len(ctitle), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (compname, len(compname), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (version, len(version), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hostname, len(hostname), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (username, len(username), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest, 1, MPI_INTEGER, 0, mpicom, ier)

    call mpi_bcast (use_lch4, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_nitrif_denitrif, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_extralakelayers, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_vichydro, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_cn, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_cndv, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_nguardrail, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_crop, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (flush_gdd20, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fertilizer, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_grainproduct, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (crop_residue_removal_frac, 1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (o3_veg_stress_method, len(o3_veg_stress_method), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (use_snicar_frc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_vancouver, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_mexicocity, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_noio, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_SSRE, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! initial file variables
    call mpi_bcast (nrevsn, len(nrevsn), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat, len(finidat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat_interp_source, len(finidat_interp_source), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat_interp_dest, len(finidat_interp_dest), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat, len(fsurdat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hillslope_file, len(hillslope_file), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmlndfrc,len(fatmlndfrc),MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (paramfile, len(paramfile) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowoptics, len(fsnowoptics),  MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowaging,  len(fsnowaging),   MPI_CHARACTER, 0, mpicom, ier)

    ! Irrigation
    call mpi_bcast(irrigate, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Crop saturated excess runoff
    call mpi_bcast(crop_fsat_equals_zero, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Whether to run tests of ncdio_pio
    call mpi_bcast(for_testing_run_ncdiopio_tests, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Various flags used for testing infrastructure for having multiple crop reproductive pools
    call mpi_bcast(for_testing_use_second_grain_pool, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(for_testing_use_repr_structure_pool, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(for_testing_no_crop_seed_replenishment, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Landunit generation
    call mpi_bcast(create_crop_landunit, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Other subgrid logic
    call mpi_bcast(run_zero_weight_urban, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(all_active, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(convert_ocean_to_land, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Number of dominant pfts and landunits. Enhance ctsm performance by
    ! reducing the number of active pfts to n_dom_pfts and
    ! active landunits to n_dom_landunits.
    ! Also choose to collapse the urban landunits to the dominant urban
    ! landunit by setting collapse_urban = .true.
    call mpi_bcast(n_dom_pfts, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast(n_dom_landunits, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast(collapse_urban, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Thresholds above which the model keeps the soil, crop, glacier, lake,
    ! wetland, and urban landunits
    call mpi_bcast(toosmall_soil, 1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast(toosmall_crop, 1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast(toosmall_glacier, 1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast(toosmall_lake, 1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast(toosmall_wetland, 1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast(toosmall_urban, 1, MPI_REAL8, 0, mpicom, ier)

    ! BGC
    call mpi_bcast (co2_type, len(co2_type), MPI_CHARACTER, 0, mpicom, ier)
    
    call mpi_bcast (use_fates, 1, MPI_LOGICAL, 0, mpicom, ier)

    if (use_cn .or. use_fates) then
       call mpi_bcast (suplnitro, len(suplnitro), MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (nfix_timeconst, 1, MPI_REAL8, 0, mpicom, ier)
       call mpi_bcast (spinup_state, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (override_bgc_restart_mismatch_dump, 1, MPI_LOGICAL, 0, mpicom, ier)
    end if

    ! isotopes
    call mpi_bcast (use_c13, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_c14, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (for_testing_allow_interp_non_ciso_to_ciso, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (fates_spitfire_mode, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (fates_harvest_mode, len(fates_harvest_mode) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_stomatal_model, len(fates_stomatal_model) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_stomatal_assimilation, len(fates_stomatal_assimilation) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_leafresp_model, len(fates_leafresp_model) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_cstarvation_model, len(fates_cstarvation_model) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_regeneration_model, len(fates_regeneration_model) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_radiation_model, len(fates_radiation_model) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_hydro_solver, len(fates_hydro_solver) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_electron_transport_model, len(fates_electron_transport_model) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (use_fates_planthydro, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_tree_damage, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_cohort_age_tracking, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_ed_st3, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_ed_prescribed_phys,  1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_inventory_init, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_fixed_biogeog, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_nocomp, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_sp, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_luh, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_lupft, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_potentialveg, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_bgc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_daylength_factor, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (fates_photosynth_acclimation, len(fates_photosynth_acclimation), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_inventory_ctrl_filename, len(fates_inventory_ctrl_filename), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_paramfile, len(fates_paramfile) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fluh_timeseries, len(fluh_timeseries) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (flandusepftdat, len(flandusepftdat) , MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (fates_parteh_mode, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (fates_seeddisp_cadence, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (fates_history_dimlevel, 2, MPI_INTEGER, 0, mpicom, ier)

    ! flexibleCN nitrogen model
    call mpi_bcast (use_flexibleCN, 1, MPI_LOGICAL, 0, mpicom, ier)
    ! TODO(bja, 2015-08) need to move some of these into a module with limited scope.
    call mpi_bcast (MM_Nuptake_opt, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (CNratio_floating, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (lnc_opt, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (reduce_dayl_factor, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (vcmax_opt, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (CN_evergreen_phenology_opt, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (carbon_resp_opt, 1, MPI_INTEGER, 0, mpicom, ier)

    call mpi_bcast (use_luna, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_soil_moisture_streams, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_excess_ice, 1, MPI_LOGICAL, 0, mpicom,ier)

    call mpi_bcast (use_lai_streams, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_cropcal_streams, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_bedrock, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_biomass_heat_storage, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_hillslope, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (downscale_hillslope_meteorology, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_hillslope_routing, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (hillslope_fsat_equals_zero, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_matrixcn, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_soil_matrixcn, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_wrt_matrixcn_diag, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (spinup_matrixcn, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (nyr_forcing, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nyr_sasu, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (iloop_avg, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (use_hydrstress, 1, MPI_LOGICAL, 0, mpicom, ier)

    if (use_cn .or. use_fates) then
       ! vertical soil mixing variables
       call mpi_bcast (som_adv_flux, 1, MPI_REAL8,  0, mpicom, ier)
       call mpi_bcast (max_depth_cryoturb, 1, MPI_REAL8,  0, mpicom, ier)

       ! C and N input vertical profiles
       call mpi_bcast (surfprof_exp,            1, MPI_REAL8,  0, mpicom, ier)
    end if

    if ((use_cn.or.use_fates) .and. use_nitrif_denitrif) then
       call mpi_bcast (no_frozen_nitrif_denitrif,  1, MPI_LOGICAL, 0, mpicom, ier)
    end if

    if (use_cn) then
       call mpi_bcast (use_c14_bombspike,  1, MPI_LOGICAL, 0, mpicom, ier)
       call mpi_bcast (atm_c14_filename,  len(atm_c14_filename), MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (use_c13_timeseries,  1, MPI_LOGICAL, 0, mpicom, ier)
       call mpi_bcast (atm_c13_filename,  len(atm_c13_filename), MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (use_fun,            1, MPI_LOGICAL, 0, mpicom, ier)
    end if

    call mpi_bcast (perchroot, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (perchroot_alt, 1, MPI_LOGICAL, 0, mpicom, ier)
    if (use_lch4) then
       call mpi_bcast (anoxia, 1, MPI_LOGICAL, 0, mpicom, ier)
    end if

    ! lakes
    call mpi_bcast (deepmixing_depthcrit,  1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (deepmixing_mixfact,    1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (lake_melt_icealb, numrad, MPI_REAL8, 0, mpicom, ier)

    ! physics variables
    call mpi_bcast (nsegspc, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (use_subgrid_fluxes , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (snow_thermal_cond_method, len(snow_thermal_cond_method), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (snow_cover_fraction_method , len(snow_cover_fraction_method), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (z0param_method , len(z0param_method), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (use_z0m_snowmelt, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (single_column,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (scmlat, 1, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (scmlon, 1, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (co2_ppmv, 1, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (albice, 2, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (soil_layerstruct_predefined,len(soil_layerstruct_predefined), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (soil_layerstruct_userdefined,size(soil_layerstruct_userdefined), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (soil_layerstruct_userdefined_nlevsoi, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (snicar_numrad_snw, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (snicar_solarspec, len(snicar_solarspec), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (snicar_dust_optics, len(snicar_dust_optics), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (snicar_snw_shape, len(snicar_snw_shape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (snicar_use_aerosol, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (snicar_snobc_intmix, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (snicar_snodst_intmix, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (do_sno_oc, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! snow pack variables
    call mpi_bcast (nlevsno, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (h2osno_max, 1, MPI_REAL8, 0, mpicom, ier)

    ! glacier_mec variables
    call mpi_bcast (maxpatch_glc, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (glc_do_dynglacier, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (glc_snow_persistence_max_days, 1, MPI_INTEGER, 0, mpicom, ier)

    ! history file variables
    call mpi_bcast (hist_empty_htapes, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, size(hist_ndens), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_avgflag_pertape, len(hist_avgflag_pertape)*size(hist_avgflag_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_type1d_pertape, max_namlen*size(hist_type1d_pertape), MPI_CHARACTER, 0, mpicom, ier)
    if (use_lch4) then
       call mpi_bcast (hist_wrtch4diag, 1, MPI_LOGICAL, 0, mpicom, ier)
    end if
    call mpi_bcast (hist_fields_list_file, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl1, max_namlen*size(hist_fexcl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl2, max_namlen*size(hist_fexcl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl3, max_namlen*size(hist_fexcl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl4, max_namlen*size(hist_fexcl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl5, max_namlen*size(hist_fexcl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl6, max_namlen*size(hist_fexcl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl7, max_namlen*size(hist_fexcl7), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl8, max_namlen*size(hist_fexcl8), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl9, max_namlen*size(hist_fexcl9), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl10,max_namlen*size(hist_fexcl10),MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl1, (max_namlen+2)*size(hist_fincl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl2, (max_namlen+2)*size(hist_fincl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl3, (max_namlen+2)*size(hist_fincl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl4, (max_namlen+2)*size(hist_fincl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl5, (max_namlen+2)*size(hist_fincl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl6, (max_namlen+2)*size(hist_fincl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl7, (max_namlen+2)*size(hist_fincl7), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl8, (max_namlen+2)*size(hist_fincl8), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl9, (max_namlen+2)*size(hist_fincl9), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl10,(max_namlen+2)*size(hist_fincl10),MPI_CHARACTER, 0, mpicom, ier)

    ! restart file variables

    call mpi_bcast (rpntfil, len(rpntfil), MPI_CHARACTER, 0, mpicom, ier)

    ! clump decomposition variables

    call mpi_bcast (clump_pproc, 1, MPI_INTEGER, 0, mpicom, ier)


  end subroutine control_spmd

  !------------------------------------------------------------------------
  subroutine control_print ()
    !
    ! !DESCRIPTION:
    ! Write out the clm namelist run control variables
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer i  !loop index
    !------------------------------------------------------------------------

    write(iulog,*) 'define run:'
    write(iulog,*) '   source                = ',trim(source)
    write(iulog,*) '   model_version         = ',trim(version)
    write(iulog,*) '   run type              = ',runtyp(nsrest+1)
    write(iulog,*) '   case title            = ',trim(ctitle)
    write(iulog,*) '   username              = ',trim(username)
    write(iulog,*) '   hostname              = ',trim(hostname)
    write(iulog,*) 'process control parameters:'
    write(iulog,*) '    use_lch4 = ', use_lch4
    write(iulog,*) '    use_nitrif_denitrif = ', use_nitrif_denitrif
    write(iulog,*) '    use_extralakelayers = ', use_extralakelayers
    write(iulog,*) '    use_vichydro = ', use_vichydro
    write(iulog,*) '    use_excess_ice = ', use_excess_ice
    write(iulog,*) '    use_cn = ', use_cn
    write(iulog,*) '    use_cndv = ', use_cndv
    write(iulog,*) '    use_crop = ', use_crop
    write(iulog,*) '    flush_gdd20 = ', flush_gdd20
    write(iulog,*) '    use_fertilizer = ', use_fertilizer
    write(iulog,*) '    use_grainproduct = ', use_grainproduct
    write(iulog,*) '    crop_residue_removal_frac = ', crop_residue_removal_frac
    write(iulog,*) '    o3_veg_stress_method = ', o3_veg_stress_method
    write(iulog,*) '    use_snicar_frc = ', use_snicar_frc
    write(iulog,*) '    snicar_use_aerosol = ',snicar_use_aerosol
    write(iulog,*) '    use_vancouver = ', use_vancouver
    write(iulog,*) '    use_mexicocity = ', use_mexicocity
    write(iulog,*) '    use_noio = ', use_noio
    write(iulog,*) '    use_SSRE = ', use_SSRE
    write(iulog,*) 'input data files:'
    write(iulog,*) '   PFT physiology and parameters file = ',trim(paramfile)
    if (fsurdat == ' ') then
       write(iulog,*) '   fsurdat, surface dataset not set'
    else
       write(iulog,*) '   surface data   = ',trim(fsurdat)
    end if
    if (hillslope_file == ' ') then
       write(iulog,*) '   hillslope_file, hillslope dataset not set'
    else
       write(iulog,*) '   hillslope data   = ',trim(hillslope_file)
    end if
    if (fatmlndfrc == ' ') then
       write(iulog,*) '   fatmlndfrc not set, setting frac/mask to 1'
    else
       write(iulog,*) '   land frac data = ',trim(fatmlndfrc)
    end if
    write(iulog,*) '   Convert ocean to land = ', convert_ocean_to_land
    write(iulog,*) '   Number of ACTIVE PFTS (0 means input pft data NOT collapsed to n_dom_pfts) =', n_dom_pfts
    write(iulog,*) '   Number of ACTIVE LANDUNITS (0 means input landunit data NOT collapsed to n_dom_landunits) =', n_dom_landunits
    write(iulog,*) '   Collapse urban landunits; done before collapsing all landunits to n_dom_landunits; .false. means do nothing i.e. keep all the urban landunits, though n_dom_landunits may still remove them =', collapse_urban
    write(iulog,*) '   Threshold above which the model keeps the soil landunit =', toosmall_soil
    write(iulog,*) '   Threshold above which the model keeps the crop landunit =', toosmall_crop
    write(iulog,*) '   Threshold above which the model keeps the glacier landunit =', toosmall_glacier
    write(iulog,*) '   Threshold above which the model keeps the lake landunit =', toosmall_lake
    write(iulog,*) '   Threshold above which the model keeps the wetland landunit =', toosmall_wetland
    write(iulog,*) '   Threshold above which the model keeps the urban landunits =', toosmall_urban

    if (use_cn .or. use_fates) then
       if (suplnitro /= suplnNon)then
          write(iulog,*) '   Supplemental Nitrogen mode is set to run over Patches: ', &
               trim(suplnitro)
       end if

       if (nfix_timeconst /= 0._r8) then
          write(iulog,*) '   nfix_timeconst, timescale for smoothing npp in N fixation term: ', nfix_timeconst
       else
          write(iulog,*) '   nfix_timeconst == zero, use standard N fixation scheme. '
       end if

       write(iulog,*) '   spinup_state, (0 = normal mode; 1 = AD spinup; 2 AAD)         : ', spinup_state
       if ( spinup_state .eq. 0 ) then
          write(iulog,*) '   model is currently NOT in AD spinup mode.'
       else if ( spinup_state .eq. 1 ) then
          write(iulog,*) '   model is currently in AD spinup mode.'
       else if ( spinup_state .eq. 2 ) then
          write(iulog,*) '   model is currently in accelerated AD spinup mode.'
       else
          call endrun(msg=' error: spinup_state can only have integer value of 0 or 1 or 2'//&
               errMsg(sourcefile, __LINE__))
       end if

       if ( use_fun ) then
          write(iulog,*) '   Fixation and Uptake of Nitrogen Model Version 2 (FUN2) is turned on for Nitrogen Competition'
       end if

       write(iulog,*) '   override_bgc_restart_mismatch_dump                     : ', override_bgc_restart_mismatch_dump
    end if

    if (use_cn .or. use_fates) then
       write(iulog, *) '   som_adv_flux, the advection term in soil mixing (m/s) : ', som_adv_flux
       write(iulog, *) '   max_depth_cryoturb (m)                                : ', max_depth_cryoturb
       write(iulog, *) '   surfprof_exp                                          : ', surfprof_exp
    end if

    if ((use_cn .or. use_fates) .and. .not. use_nitrif_denitrif) then
       write(iulog, *) '   no_frozen_nitrif_denitrif                             : ', no_frozen_nitrif_denitrif
    end if

    if (use_cn) then
       write(iulog, *) '  use_c13                                                : ', use_c13
       write(iulog, *) '  use_c13_timeseries                                     : ', use_c13_timeseries
       write(iulog, *) '  atm_c13_filename                                       : ', atm_c13_filename
       write(iulog, *) '  use_c14                                                : ', use_c14
       write(iulog, *) '  use_c14_bombspike                                      : ', use_c14_bombspike
       write(iulog, *) '  atm_c14_filename                                       : ', atm_c14_filename
       write(iulog, *) '  for_testing_allow_interp_non_ciso_to_ciso              : ', for_testing_allow_interp_non_ciso_to_ciso
    end if

    if (fsnowoptics == ' ') then
       write(iulog,'(a)') '   snow optical properties file NOT set'
    else
       write(iulog,'(a)') '   snow optical properties file = '//trim(fsnowoptics)
    endif
    if (fsnowaging == ' ') then
       write(iulog,'(a)') '   snow aging parameters file NOT set'
    else
       write(iulog,'(a)') '   snow aging parameters file = '//trim(fsnowaging)
    endif

    write(iulog,*) '   SNICAR: downward solar radiation spectrum type =', snicar_solarspec
    write(iulog,*) '   SNICAR: dust optics type = ', snicar_dust_optics
    write(iulog,*) '   SNICAR: number of bands in snow albedo calculation =', snicar_numrad_snw
    write(iulog,*) '   SNICAR: snow grain shape type = ', snicar_snw_shape
    write(iulog,*) '   SNICAR: BC-snow internal mixing = ', snicar_snobc_intmix
    write(iulog,*) '   SNICAR: dust-snow internal mixing = ', snicar_snodst_intmix
    write(iulog,*) '   SNICAR: OC in snow = ', do_sno_oc

    write(iulog,'(a,i8)') '   Number of snow layers =', nlevsno
    write(iulog,'(a,d20.10)') '   Max snow depth (mm) =', h2osno_max

    write(iulog,'(a,i8)') '   glc number of elevation classes =', maxpatch_glc

    if (glc_do_dynglacier) then
       write(iulog,*) '   glc CLM glacier areas and topography WILL evolve dynamically'
    else
       write(iulog,'(a)') '   glc CLM glacier areas and topography will NOT evolve dynamically'
    end if
    write(iulog,'(a,i8)') '   glc snow persistence max days = ', glc_snow_persistence_max_days

    if (nsrest == nsrStartup) then
       if (finidat /= ' ') then
          write(iulog,'(a)') '   initial data: ', trim(finidat)
       else if (finidat_interp_source /= ' ') then
          write(iulog,'(a)') '   initial data interpolated from: '// trim(finidat_interp_source)
       else
          write(iulog,'(a)') '   initial data created by model (cold start)'
       end if
    else
       write(iulog,*) '   restart data   = ',trim(nrevsn)
    end if

    write(iulog,'(a)') '   atmospheric forcing data is from cesm atm model'
    write(iulog,'(a)') 'Restart parameters:'
    write(iulog,'(a)')'   restart pointer file directory     = '//trim(rpntdir)
    write(iulog,'(a)')'   restart pointer file name          = '//trim(rpntfil)
    write(iulog,'(a)') 'model physics parameters:'

    if ( trim(co2_type) == 'constant' )then
       write(iulog,'(a,d20.10)') '   CO2 volume mixing ratio   (umol/mol)   = ', co2_ppmv
    else
       write(iulog,'(a)'       ) '   CO2 volume mixing ratio                = '//trim(co2_type)
    end if

    write(iulog,*) '   land-ice albedos      (unitless 0-1)   = ', albice
    write(iulog,*) '   hillslope hydrology    = ', use_hillslope
    write(iulog,*) '   downscale hillslope meteorology  = ', downscale_hillslope_meteorology
    write(iulog,*) '   hillslope routing      = ', use_hillslope_routing
    write(iulog,*) '   hillslope_fsat_equals_zero       = ', hillslope_fsat_equals_zero
    write(iulog,*) '   pre-defined soil layer structure = ', soil_layerstruct_predefined
    write(iulog,*) '   user-defined soil layer structure = ', soil_layerstruct_userdefined
    write(iulog,*) '   user-defined number of soil layers = ', soil_layerstruct_userdefined_nlevsoi
    write(iulog,*) '   plant hydraulic stress = ', use_hydrstress
    if (nsrest == nsrContinue) then
       write(iulog,*) 'restart warning:'
       write(iulog,*) '   Namelist not checked for agreement with initial run.'
       write(iulog,*) '   Namelist should not differ except for ending time step and run type'
    end if
    if (nsrest == nsrBranch) then
       write(iulog,'(a)') 'branch warning:'
       write(iulog,'(a)') '   Namelist not checked for agreement with initial run.'
       write(iulog,'(a)') '   Surface data set and reference date should not differ from initial run'
    end if
    write(iulog,*) '   nsegspc              = ',nsegspc
    ! New fields
    write(iulog,*) ' perchroot (plant water stress based on unfrozen layers only) = ',perchroot
    write(iulog,*) ' perchroot (plant water stress based on time-integrated active layer only) = ',perchroot
    if (use_lch4) then
       write(iulog,*) ' anoxia (applied to soil decomposition)             = ',anoxia
    end if
    ! Lakes
    write(iulog,*)
    write(iulog,*) 'Lake Model Namelists:'
    write(iulog,*) 'Increased mixing relative to Hostetler wind-driven eddy expression ',&
                   'will be used for deep lakes exceeding depth ', deepmixing_depthcrit,&
                      ' by a factor of ', deepmixing_mixfact, '.'
    write(iulog,*) 'Albedo over melting lakes will approach values (visible, NIR):', lake_melt_icealb, &
                   'as compared with 0.60, 0.40 for cold frozen lakes with no snow.'

    write(iulog, *) 'plant nitrogen model namelists:'
    write(iulog, *) '  use_flexibleCN = ', use_flexibleCN
    if (use_flexibleCN) then
       write(iulog, *) '    MM_Nuptake_opt = ', MM_Nuptake_opt
       write(iulog, *) '    CNratio_floating = ', CNratio_floating
       write(iulog, *) '    lnc_opt = ', lnc_opt
       write(iulog, *) '    reduce_dayl_factor = ', reduce_dayl_factor
       write(iulog, *) '    vcmax_opt = ', vcmax_opt
       write(iulog, *) '    CN_evergreen_phenology_opt = ', CN_evergreen_phenology_opt
       write(iulog, *) '    carbon_resp_opt = ', carbon_resp_opt
    end if
    write(iulog, *) '  use_luna = ', use_luna

    write(iulog, *) '  ED/FATES: '
    write(iulog, *) '    use_fates = ', use_fates
    if (use_fates) then
       write(iulog, *) '    fates_spitfire_mode = ', fates_spitfire_mode
       write(iulog, *) '    fates_harvest_mode = ', fates_harvest_mode
       write(iulog, *) '    fates_stomatal_model = ', fates_stomatal_model
       write(iulog, *) '    fates_stomatal_assimilation = ', fates_stomatal_assimilation
       write(iulog, *) '    fates_leafresp_model = ', fates_leafresp_model
       write(iulog, *) '    fates_cstarvation_model = ', fates_cstarvation_model
       write(iulog, *) '    fates_regeneration_model = ', fates_regeneration_model
       write(iulog, *) '    fates_radiation_model = ', fates_radiation_model
       write(iulog, *) '    fates_hydro_solver = ', fates_hydro_solver
       write(iulog, *) '    fates_electron_transport_model = ', fates_electron_transport_model
       write(iulog, *) '    fates_paramfile = ', fates_paramfile
       write(iulog, *) '    fates_parteh_mode = ', fates_parteh_mode
       write(iulog, *) '    fates_photosynth_acclimation = ', trim(fates_photosynth_acclimation)
       write(iulog, *) '    use_fates_planthydro = ', use_fates_planthydro
       write(iulog, *) '    use_fates_tree_damage = ', use_fates_tree_damage
       write(iulog, *) '    use_fates_cohort_age_tracking = ', use_fates_cohort_age_tracking
       write(iulog, *) '    use_fates_daylength_factor = ', use_fates_daylength_factor
       write(iulog, *) '    use_fates_ed_st3 = ',use_fates_ed_st3
       write(iulog, *) '    use_fates_ed_prescribed_phys = ',use_fates_ed_prescribed_phys
       write(iulog, *) '    use_fates_inventory_init = ',use_fates_inventory_init
       write(iulog, *) '    use_fates_fixed_biogeog = ', use_fates_fixed_biogeog
       write(iulog, *) '    use_fates_nocomp = ', use_fates_nocomp
       write(iulog, *) '    use_fates_sp = ', use_fates_sp
       write(iulog, *) '    use_fates_luh= ', use_fates_luh
       write(iulog, *) '    use_fates_lupft= ', use_fates_lupft
       write(iulog, *) '    use_fates_potentialveg = ', use_fates_potentialveg
       write(iulog, *) '    fluh_timeseries = ', trim(fluh_timeseries)
       write(iulog, *) '    flandusepftdat = ', trim(flandusepftdat)
       write(iulog, *) '    fates_seeddisp_cadence = ', fates_seeddisp_cadence
       write(iulog, *) '    fates_seeddisp_cadence: 0, 1, 2, 3 => off, daily, monthly, or yearly dispersal'
       write(iulog, *) '    fates_inventory_ctrl_filename = ', trim(fates_inventory_ctrl_filename)
    end if
  end subroutine control_print


  !-----------------------------------------------------------------------
  subroutine check_missing_initdata_status(finidat_interp_dest)
   !
   ! !DESCRIPTION:
   ! Checks that the finidat_interp_dest .status file was written (i.e., that write of
   ! finidat_interp_dest succeeded)
   !
   ! !ARGUMENTS:
   character(len=*), intent(in)    :: finidat_interp_dest
   !
   ! !LOCAL VARIABLES:
   logical                    :: lexists
   integer                    :: klen
   character(len=fname_len)   :: status_file
   character(len=*), parameter :: subname = 'check_missing_initdata_status'
   !-----------------------------------------------------------------------

    klen = len_trim(finidat_interp_dest) - 3 ! remove the .nc
    status_file = finidat_interp_dest(1:klen)//'.status'
    inquire(file=trim(status_file), exist=lexists)
    if (.not. lexists) then
       if (masterproc) then
          write(iulog,'(a)')' failed to find file '//trim(status_file)
          write(iulog,'(a)')' this indicates a problem in creating '//trim(finidat_interp_dest)
          write(iulog,'(a)')' remove '//trim(finidat_interp_dest)//' and try again'
       end if
       call endrun(subname//': finidat_interp_dest file exists but is probably bad')
    end if

  end subroutine check_missing_initdata_status


  !-----------------------------------------------------------------------
  subroutine apply_use_init_interp(finidat_interp_dest, finidat, finidat_interp_source)
    !
    ! !DESCRIPTION:
    ! Applies the use_init_interp option, setting finidat_interp_source to finidat
    !
    ! Should be called if use_init_interp is true.
    !
    ! Does error checking to ensure that it is valid to set use_init_interp to true,
    ! given the values of finidat and finidat_interp_source.
    !
    ! !USES:
    use netcdf
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)    :: finidat_interp_dest
    character(len=*), intent(inout) :: finidat
    character(len=*), intent(inout) :: finidat_interp_source
    !
    ! !LOCAL VARIABLES:
    logical                    :: lexists
    integer                    :: ncid
    character(len=fname_len)   :: initial_source_file
    integer                    :: status
    character(len=*), parameter :: subname = 'apply_use_init_interp'
    !-----------------------------------------------------------------------

    if (finidat == ' ') then
       write(iulog,*)' WARNING: Setting use_init_interp has no effect if finidat is not also set'
    end if

    if (finidat_interp_source /= ' ') then
       call endrun(msg=' ERROR: Cannot set use_init_interp if finidat_interp_source is &
            &already set')
    end if

    if ( get_filename(finidat) == &
         get_filename(finidat_interp_dest)) then
       write(iulog,*) 'ERROR: With use_init_interp, cannot use the same filename for source and dest'
       write(iulog,*) '(because this will lead to the source being overwritten before it is read).'
       write(iulog,*) 'finidat             = ', trim(get_filename(finidat))
       write(iulog,*) 'finidat_interp_dest = ', trim(get_filename(finidat_interp_dest))
       write(iulog,*) '(Even if the two files are in different directories, you cannot use the same filename'
       write(iulog,*) 'due to problems related to <https://github.com/ESCOMP/ctsm/issues/329>.)'
       write(iulog,*) 'As a workaround, copy or move the finidat file to a different name.'
       call endrun(msg=' ERROR: With use_init_interp, cannot use the same filename for source and dest')
    end if

    inquire(file=trim(finidat_interp_dest), exist=lexists)
    if (lexists) then

       ! Check that the status file also exists (i.e., that finidat_interp_dest was written successfully)
       call check_missing_initdata_status(finidat_interp_dest)

       ! open the input file and check for the name of the input source file
       status = nf90_open(trim(finidat_interp_dest), 0, ncid)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_get_att(ncid, NF90_GLOBAL, 'initial_source_file', initial_source_file)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_close(ncid)
       if (status /= nf90_noerr) call handle_err(status)

       ! If the input source file in finidat_interp_dest does not match the namelist input for mapping
       ! then use the namelist input (i.e. finidat) - otherwise just use the generated file
       if (trim(initial_source_file) /= trim(finidat)) then
          finidat_interp_source = finidat
          finidat = ' '
          write(iulog,'(a)')' interpolating initial dataset using source file '//trim(finidat_interp_source)
       else
          finidat = trim(finidat_interp_dest)
          write(iulog,'(a)')' using interpolated initial dataset file '//trim(finidat)
       end if
    else
       finidat_interp_source = finidat
       finidat = ' '
       write(iulog,'(a)')' interpolating initial dataset using source file '//trim(finidat_interp_source)
    end if

  contains
    subroutine handle_err(status)
      integer, intent ( in) :: status
      if (status /= nf90_noerr) then
         if (masterproc) then
            write(iulog,'(a)') trim(nf90_strerror(status))
         end if
         call endrun()
      end if
    end subroutine handle_err

  end subroutine apply_use_init_interp


end module controlMod
