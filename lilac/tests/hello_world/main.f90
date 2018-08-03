program main

  ! modules
  use ESMF
  ! use lilac, ONLY : lilac_init

  implicit none

  ! local variables
  integer:: rc

  ! call lilac_init()
  ! TODO fix linking with lilac
  call ESMF_Initialize(rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  print *, "Hello LILAC World"

  call ESMF_Finalize()

end program main
