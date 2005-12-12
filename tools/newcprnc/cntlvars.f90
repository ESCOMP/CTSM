module cntlvars
!
! Control variables set in main program and used in lower level routines.
! All are available as cmd-line args so set default values here.
!
! Indices for printing. Init to zero so can test on > 0 for valid user input
!
   integer :: iprs = 0            ! start index 1st dimension
   integer :: ipre = 0            ! end index 1st dimension
   integer :: jprs = 0            ! start index 2nd dimension (if applicable)
   integer :: jpre = 0            ! end index 2nd dimension (if applicable)
   integer :: kprs = 0            ! start index 3rd dimension (if applicable)
   integer :: kpre = 0            ! end index 3rd dimension (if applicable)
!
! Logical variables
!
   logical :: verbose = .false.   ! verbose output (currently unused)
   logical :: lnorm   = .false.   ! print level by level L2 and L-INF norm stats (currently unused)
   logical :: constps = .false.   ! always use constant surface pressure (currently unused)
   logical :: matchts = .true.    ! Align timesteps before comparison (default true)
end module cntlvars
