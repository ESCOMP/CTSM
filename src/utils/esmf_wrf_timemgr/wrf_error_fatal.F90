subroutine wrf_error_fatal(msg)
  use abortutils
  implicit none
  character(len=*), intent(in) :: msg
  call endrun(msg)
end subroutine wrf_error_fatal
