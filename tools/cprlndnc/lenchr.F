integer function lenchr(chrstg)
  implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
  character*(*) chrstg       !  Input character string
!
! Local workspace
!
  integer l
!
!-----------------------------------------------------------------------
!
  lenchr = 0
  do l=len(chrstg),1,-1
    if (chrstg(l:l).ne.' ' .and. chrstg(l:l).ne.char(0)) then
      lenchr = l
      goto 10
    end if
  end do
10 return
end function lenchr


