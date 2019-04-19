program demo_lilac_driver

  ! modules
  use ESMF
  use LilacMod

  implicit none

  real, dimension(10) :: dum_var1
  real, dimension(10) :: dum_var2

  !real, dimension(10) :: t_phy   ! temperature (K)
  !real, dimension(10) :: th_phy  ! potential temperature (K)
  !real, dimension(10,10) :: rho     !

  call random_number(dum_var1)
  call random_number(dum_var2)

  !call random_number(t_phy)
  !call random_number(th_phy)

  print *, "dum_var1 = ", dum_var1
  print *, "dum_var2 = ", dum_var2
 
  call lilac_init(dum_var1, dum_var2)
  call lilac_run(dum_var1, dum_var2)

  call ESMF_Finalize()


end program demo_lilac_driver

