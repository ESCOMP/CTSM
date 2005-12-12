integer function ismax(n,sx,incx)
  real*8 sx(*)
  real*8 max
  integer imax
  max = sx(1)
  imax = 1
  do i=1+incx,n,incx
    if (max .lt. sx(i)) then
      max = sx(i)
      imax = i
    endif
  end do
  ismax = imax
  return
end function ismax
