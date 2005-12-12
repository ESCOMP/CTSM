subroutine prhddiff (nlev, nlon, nlat)

  use precision
  use header
  implicit none

  ! Input arguments

  integer, intent(in) :: nlev, nlon, nlat

  ! Local workspace

  integer k, i, j, n
  !---------------------------------------------------------------

  write(6,*)'SUMMARY OF IMPORTANT HEADER DIFFERENCES:'

  do n=1,ntime(1)
     if (nsteph(n,1) /= nsteph(n,2)) write(6,100) 'nsteph:', nsteph(n,1), nsteph(n,2)
     if (ncdate(n,1) /= ncdate(n,2)) write(6,100) 'ncdate:', ncdate(n,1), ncdate(n,2)
     if (ncsec(n,1)  /= ncsec(n,2) ) write(6,100) 'ncsec :', ncsec(n,1) , ncsec(n,2) 
  end do

  write(6,*)
  do j=1,nlat
     do i=1,nlon
        if (area(i,j,1) /= area(i,j,2)) then
           write(6,200) 'area(',i,j,'):',area(i,j,1), area(i,j,2)
           stop
        endif
     end do
  end do

  return
100 format(a,2i9)
200 format(a,i3,a,1p,2e23.15)

end subroutine prhddiff
