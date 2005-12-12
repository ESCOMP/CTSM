module specialvalues
!
! Special values used in various initializations
!
   use prec, only: r8

   integer, parameter :: uninit_int = -99999
   
   real(r8), parameter :: inf = O'0777600000000000000000'  ! infinity (assumes IEEE)
   real(r8), parameter :: uninit_r8 = inf
end module specialvalues
