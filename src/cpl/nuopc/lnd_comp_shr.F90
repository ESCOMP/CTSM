module lnd_comp_shr

  ! Model mesh info is here in order to be leveraged by CDEPS in line calls

  use ESMF        , only : ESMF_Clock, ESMF_Mesh
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl  , only : fl => fname_len

  implicit none
  public

  type(ESMF_Clock)  :: model_clock    ! model clock
  type(ESMF_Mesh)   :: mesh           ! model_mesh
  character(len=fl) :: model_meshfile ! model mesh file

end module lnd_comp_shr
