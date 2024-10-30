module CNSharedParamsMod

  !-----------------------------------------------------------------------
  !
  ! Parameters that are shared by the Carbon Nitrogen Biogeochemistry modules
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none

  ! CNParamsShareInst.  PGI wants the type decl. public but the instance
  ! is indeed protected.  A generic private statement at the start of the module
  ! overrides the protected functionality with PGI

  type, public  :: CNParamsShareType
      real(r8) :: Q10                   ! temperature dependence
      real(r8) :: minpsi                ! minimum soil water potential for heterotrophic resp	  
      real(r8) :: maxpsi                ! maximum soil water potential for heterotrophic resp
      real(r8) :: rf_cwdl2              ! respiration fraction in CWD to litter2 transition (frac)
      real(r8) :: tau_cwd               ! corrected fragmentation rate constant CWD, century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1 (1/0.3) (1/yr)
      real(r8) :: cwd_flig              ! lignin fraction of coarse woody debris
      real(r8) :: froz_q10              ! separate q10 for frozen soil respiration rates
      real(r8) :: decomp_depth_efolding ! e-folding depth for reduction in decomposition (m) 
      real(r8) :: mino2lim              ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
      real(r8) :: organic_max           ! organic matter content (kg/m3) where soil is assumed to act like peat
      logical  :: constrain_stress_deciduous_onset ! if true use additional constraint on stress deciduous onset trigger
  end type CNParamsShareType

  type(CNParamsShareType), protected :: CNParamsShareInst

  ! Public subroutines
  public :: CNParamsReadShared              ! Read in CN shared parameters
  public :: CNParamsSetSoilDepth            ! Set the soil depth needed for CNPhenology
  public :: CNParamsReadShared_namelist     ! Read in CN shared namelist items

  ! Public data

  logical, public :: use_matrixcn = .false.  ! true => use cn matrix solution
  logical, public :: use_fun      = .false.             ! Use the FUN2.0 model
  integer, public :: nlev_soildecomp_standard = 5
  integer, public :: upper_soil_layer = -1              ! Upper soil layer to use for 10-day average in CNPhenology

  ! Private subroutines and data

  private :: CNParamsReadShared_netcdf                  ! Read shared parameters from NetCDF file

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  subroutine CNParamsReadShared(ncid, namelist_file)

    use ncdio_pio      , only : file_desc_t

    type(file_desc_t), intent(inout) :: ncid   ! pio netCDF file id
    character(len=*),  intent(in) :: namelist_file

    call CNParamsReadShared_netcdf(ncid)
    call CNParamsReadShared_namelist(namelist_file)

  end subroutine CNParamsReadShared
  
  !-----------------------------------------------------------------------

  subroutine CNParamsSetSoilDepth( )
    use initVerticalMod, only : find_soil_layer_containing_depth
    ! Set the soil depth needed for CNPhenology
    call find_soil_layer_containing_depth ( 0.12_r8, upper_soil_layer )
  end subroutine CNParamsSetSoilDepth
  !-----------------------------------------------------------------------
  subroutine CNParamsReadShared_netcdf(ncid)
    !
    use ncdio_pio   , only : file_desc_t, ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'CNParamsReadShared'
    character(len=100) :: errCode = '-Error reading in CN and BGC shared params file. Var:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! netcdf read here
    !
    tString='q10_mr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%Q10=tempr

    tString='minpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%minpsi=tempr 

    tString='maxpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%maxpsi=tempr

    tString='rf_cwdl2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%rf_cwdl2=tempr

    tString='tau_cwd'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%tau_cwd=tempr

    tString='cwd_flig'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%cwd_flig=tempr 

    tString='decomp_depth_efolding'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%decomp_depth_efolding=tempr

    tString='froz_q10'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%froz_q10=tempr   

    tString='mino2lim'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%mino2lim=tempr 
    !CNParamsShareInst%mino2lim=0.2_r8 

    tString='organic_max'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    CNParamsShareInst%organic_max=tempr

  end subroutine CNParamsReadShared_netcdf
  
  !-----------------------------------------------------------------------
  subroutine CNParamsReadShared_namelist(namelist_file)
    !
    ! !DESCRIPTION:
    ! Read and initialize CN Shared parameteres from the namelist.
    !
    ! !USES:
    use fileutils   , only : relavu, getavu
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_LOGICAL
    use shr_nl_mod  , only : shr_nl_find_group_name
    use shr_log_mod , only : errMsg => shr_log_errMsg
    use clm_varctl  , only : iulog
    use abortutils  , only : endrun
    use shr_mpi_mod , only : shr_mpi_bcast
    
    !
    implicit none
    !

    character(len=*), intent(in) :: namelist_file
    
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    logical  :: constrain_stress_deciduous_onset = .false.

    character(len=32) :: subroutine_name = 'CNParamsReadNamelist'
    character(len=10) :: namelist_group = 'bgc_shared'

    !-----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    namelist /bgc_shared/ &
         constrain_stress_deciduous_onset


    ! Read namelist from standard input.
    if (masterproc) then

       write(iulog,*) 'Attempting to read CN/BGC shared namelist parameters .....'
       unitn = getavu()
       write(iulog,*) 'Read in ' // namelist_group // ' namelist from: ', trim(namelist_file)
       open( unitn, file=trim(namelist_file), status='old' )
       call shr_nl_find_group_name(unitn, namelist_group, status=ierr)
       if (ierr == 0) then
          read(unitn, bgc_shared, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='error in reading in ' // namelist_group // ' namelist' // &
                  errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='error in finding ' // namelist_group // ' namelist' // &
                  errMsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )

    end if ! masterproc

    ! Broadcast the parameters from master
    call shr_mpi_bcast ( constrain_stress_deciduous_onset, mpicom )

    ! Save the parameter to the instance
    CNParamsShareInst%constrain_stress_deciduous_onset = constrain_stress_deciduous_onset

    ! Output read parameters to the lnd.log
    if (masterproc) then
       write(iulog,*) 'CN/BGC shared namelist parameters:'
       write(iulog,*)' '
       write(iulog,*)'  constrain_stress_deciduous_onset = ',constrain_stress_deciduous_onset

       write(iulog,*)

    end if

  end subroutine CNParamsReadShared_namelist

end module CNSharedParamsMod
