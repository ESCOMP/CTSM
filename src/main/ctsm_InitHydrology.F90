subroutine ctsm_InitHydrology( NLFilename )
!
!DESCRIPTION
! Initialize implementation methods for different hydrology sub-modules 
! This is created for unit-based sensitivity tests
! created by Jinyun Tang, Mar 22, 2014.

  ! !USES:
  use ctsm_Spmd       , only : masterproc, mpicom
  use ctsm_FileUtils     , only : getavu, relavu, opnfil
  use shr_nl_mod    , only : shr_nl_find_group_name
  use shr_mpi_mod   , only : shr_mpi_bcast
    
  use ctsm_FuncPedotransfer,  only : init_pedof
  use ctsm_RootBiophys,       only : init_rootprof
  use ctsm_SoilWaterMovement, only : init_soilwater_movement
  use ctsm_SoilMoistStress,   only : init_root_moist_stress
implicit none

  character(len=*), intent(IN) :: NLFilename ! Namelist filename

  !In future versions, a namelist will be created here to
  !set up options for different sub-models, the namelist file
  !will also be passed into this different initializing methods
  !to read in their local parameters, Jinyun Tang, Mar 29, 2014
  
  call init_pedof

  call init_rootprof(NLFilename)
  
  call init_soilwater_movement

  call init_root_moist_stress
  
! remove due to circular dependency of nlfilename, read namlist 
! in controlmod instead, as is done for canopyhydrology
!  call init_soil_resistance
  
end subroutine ctsm_InitHydrology
