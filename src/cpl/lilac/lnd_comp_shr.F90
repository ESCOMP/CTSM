module lnd_comp_shr

  ! Model mesh info is here in order to be leveraged by CDEPS in line calls

  use ESMF        , only : ESMF_Clock, ESMF_Mesh
  use shr_kind_mod, only : r8 => shr_kind_r8, cl=>shr_kind_cl

  implicit none
  public

  type(ESMF_Clock)  :: model_clock    ! model clock
  type(ESMF_Mesh)   :: mesh           ! model_mesh
  character(len=cl) :: model_meshfile ! model mesh file 

end module lnd_comp_shr
