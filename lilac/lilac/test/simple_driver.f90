program simple_driver
  implicit none
  integer :: t

  t = 1

  if (t == 1) then
     write(*,*) "on this line"
  else
     write(*,*) "but not here"
  end if
end program
