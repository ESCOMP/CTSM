module ch4varcon

  !-----------------------------------------------------------------------
  ! Module containing CH4 parameters and logical switches and routine to read constants from CLM namelist.
  !
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use clm_varctl  , only : NLFileName_in
  !
  implicit none
  !
  ! Methane Model Parameters
  !
  private

  logical, public :: use_aereoxid_prog = .true. ! if false then aereoxid is read off of
  ! the parameter file and may be modifed by the user (default aereoxid on the
  ! file is 0.0).

  logical, public :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                      ! Transpiration loss assumes that the methane concentration in dissolved soil
                                      ! water remains constant through the plant and is released when the water evaporates
                                      ! from the stomata.
                                      ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

  logical, public :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                     ! (Methane can be produced, and CO2 produced from methane oxidation,
                                     ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                     ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                     ! will still be averaged into the concentration _sat history fields.

  logical, public :: usephfact = .false. ! Switch to use pH factor in methane production

  logical, public :: replenishlakec = .true. ! Switch for keeping carbon storage under lakes constant
                                      ! so that lakes do not affect the carbon balance
                                      ! Good for long term rather than transient warming experiments
               ! NOTE SWITCHING THIS OFF ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
               ! IN THIS MODE.

  ! inundatrion fraction -- which is used in methane code and potentially soil code
  integer, public :: finundation_mtd        ! Finundation method type to use, one of the following
  integer, public, parameter :: finundation_mtd_h2osfc        = 0 ! Use prognostic fsat h2osfc
  integer, public, parameter :: finundation_mtd_ZWT_inversion = 1 ! Use inversion of ZWT to Prigent satellite inundation obs. data
  integer, public, parameter :: finundation_mtd_TWS_inversion = 2 ! Use inversion of TWS to Prigent satellite inundation obs. data

  logical, public :: usefrootc = .false.    ! Use CLMCN fine root C rather than ann NPP & LAI based parameterization to
                                    ! calculate tiller C for aerenchyma area calculation.
                                    ! The NPP & LAI param. was based on Wania for Arctic sedges and may not be
                                    ! appropriate for woody Patches, although nongrassporosratio above partly adjusts
                                    ! for this.  However, using fine root C reduces the aerenchyma area by a large
                                    ! factor.

  logical, public :: ch4offline = .true.    ! true --> Methane is not passed between the land & atmosphere.
                                    ! NEM is not added to NEE flux to atm. to correct for methane production,
                                    ! and ambient CH4 is set to constant 2009 value.

  logical, public :: ch4rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                    ! methanogenesis.
                                    ! Note: this option has not been extensively tested.
                                    ! Currently hardwired off.

  logical, public :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                        ! Currently hardwired off; expression is crude.

  logical, public :: ch4frzout = .false.    ! Exclude CH4 from frozen fraction of soil pore H2O, to simulate "freeze-out" pulse
                                    ! as in Mastepanov 2008.
                                    ! Causes slight increase in emissions in the fall and decrease in the spring.
                                    ! Currently hardwired off; small impact.

  public :: ch4conrd ! Read and initialize CH4 constants

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ch4conrd ()
    !
    ! !DESCRIPTION:
    ! Read and initialize CH4 constants
    !
    ! !USES:
    use fileutils   , only : relavu, getavu
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_LOGICAL, MPI_INTEGER
    use shr_nl_mod  , only : shr_nl_find_group_name
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    implicit none
    !
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'ch4conrd'  ! subroutine name
    character(len=50) :: finundation_method = 'ZWT_inversion'    ! String for finundation_method to use
    !-----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! Driver
    namelist /ch4par_in/ &
         ch4offline, replenishlakec, allowlakeprod, finundation_method

    ! Production
    namelist /ch4par_in/ &
         usephfact 

    ! Methane
    namelist /ch4par_in/ &
         use_aereoxid_prog, usefrootc

       ! ----------------------------------------------------------------------
       ! Read namelist from standard input.
       ! ----------------------------------------------------------------------

    if (masterproc) then

       write(iulog,*) 'Attempting to read CH4 parameters .....'
       unitn = getavu()
       write(iulog,*) 'Read in ch4par_in namelist from: ', trim(NLFilename_in)
       open( unitn, file=trim(NLFilename_in), status='old' )
       call shr_nl_find_group_name(unitn, 'ch4par_in', status=ierr)
       if (ierr == 0) then
          read(unitn, ch4par_in, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='error in reading in ch4par_in namelist'//&
                  errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='error in finding ch4par_in namelist'//&
                  errMsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )

       if (      trim(finundation_method) == "h2osfc"        )then
         finundation_mtd = finundation_mtd_h2osfc
       else if ( trim(finundation_method) == "TWS_inversion" )then
         finundation_mtd = finundation_mtd_tws_inversion
       else if ( trim(finundation_method) == "ZWT_inversion" )then
         finundation_mtd = finundation_mtd_zwt_inversion
       else
          call endrun(msg='error bad value for finundation_method in ch4par_in namelist'//&
                  errMsg(sourcefile, __LINE__))
       end if

    end if ! masterproc


    call mpi_bcast ( use_aereoxid_prog, 1 , MPI_LOGICAL, 0, mpicom, ierr )          
    call mpi_bcast (allowlakeprod,      1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (usephfact,          1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (replenishlakec,     1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (finundation_mtd,    1 , MPI_INTEGER, 0, mpicom, ierr)            
    call mpi_bcast (usefrootc,          1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (ch4offline,         1 , MPI_LOGICAL, 0, mpicom, ierr)            

    if (masterproc) then
       write(iulog,*) 'Successfully read CH4 namelist'
       write(iulog,*)' '
       write(iulog,*)'allowlakeprod      = ', allowlakeprod
       write(iulog,*)'usephfact          = ', usephfact
       write(iulog,*)'replenishlakec     = ', replenishlakec
       write(iulog,*)'finundation_method = ', finundation_method
       write(iulog,*)'usefrootc          = ', usefrootc
       write(iulog,*)'ch4offline         = ', ch4offline
       write(iulog,*)'use_aereoxid_prog  = ', use_aereoxid_prog 

       if (ch4offline) write(iulog,*)'CH4 Model will be running offline and not affect fluxes to atmosphere'

       if ( .not. use_aereoxid_prog ) then
          write(iulog,*) 'Aerenchyma oxidation (aereoxid) value is being read from '//&
            ' the parameters file'
       endif

       if (.not. allowlakeprod) write(iulog,*) 'Lake production has been disabled. '// &
          '  Lakes will not factor into CH4 BGC.  "Sat" history fields will not average'// &
          '  over lakes except for concentrations, which will average zero from lakes.'
       if (.not. replenishlakec .and. .not. ch4offline) write(iulog,*)'LAKE SOIL CARBON WILL NOT BE REPLENISHED BUT INSTEAD ',&
             'WILL BE TRANSIENTLY RELEASED: COUPLED MODEL WILL NOT CONSERVE CARBON IN THIS MODE!'
       write(iulog,*)'Successfully initialized CH4 parameters from namelist.'
       write(iulog,*)

    end if

  end subroutine ch4conrd

end module ch4varcon

