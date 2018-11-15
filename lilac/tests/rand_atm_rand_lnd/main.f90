program main

  use rand_test, only : atm_driver
  use ESMF

  implicit none

  ! local variables
  integer:: rc

  rc = 0
  print *, "Running Atmosphere Driver"

  call atm_driver(rc)

  if (rc /= ESMF_SUCCESS) stop 1

  print *, "Done Running Atmosphere Driver"

end program main
