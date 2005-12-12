subroutine printstats (twocases, name, numlev, nlev, &
                       nlon, nlat, wsum)
  use precision
  use stats

  implicit none
  !
  ! Arguments
  !
  logical twocases
  character*(*), intent(in) :: name
  integer, intent(in) :: numlev, nlev, nlon, nlat
  real(r8), intent(in) :: wsum
  !
  ! Local workspace
  !
  real(r8) digworst
  real(r8) digbar
  !
  ! Externals
  !
  integer lenchr
  external lenchr

  if (twocases) then  ! Two cases

     if (ndif.gt.0) then    ! Print difference info

        if (numlev.gt.1) then
           write(6,810)imx(1),jmx(1),kmx(1), &
                imn(1),jmn(1),kmn(1), &
                idmxsv ,jdmxsv ,kdmxsv , &
                irdmxsv,jrdmxsv,krdmxsv
           write(6,803)name,ndif,npos,xmx(1),xmn(1), &
                difmx,dmxsv(1),rdifmx,rdmxsv(1), &
                xmx(2),xmn(2),dmxsv(2),rdmxsv(2)
           write(6,812)imx(2),jmx(2),kmx(2), &
                imn(2),jmn(2),kmn(2)
        else
           write(6,811)imx(1),jmx(1), &
                imn(1),jmn(1), &
                idmxsv ,jdmxsv , &
                irdmxsv,jrdmxsv
           write(6,803)name,ndif,npos,xmx(1),xmn(1), &
                difmx,dmxsv(1),rdifmx,rdmxsv(1), &
                xmx(2),xmn(2),dmxsv(2),rdmxsv(2)
           write(6,813)imx(2),jmx(2), &
                imn(2),jmn(2)
        end if

        if (npos > 0) then
           rms = sqrt(rms/npos)
           xbar(1) = xbar(1)/npos
           xbar(2) = xbar(2)/npos
           rdbar = rdbar/npos
        end if

        ! Compute # digits accuracy for worst case & avg differences

        digbar = rdlnbar/ndif
        digworst = log10(1./rdifmx)

        write(6,805) xbar(1),rms,rdbar, &
             xbar(2),digbar,digworst

     else           ! Zero differences

        if (npos > 0) then
           xbar(1) = xbar(1)/npos
           xbar(2) = xbar(2)/npos
        end if

        if (numlev.gt.1) then
           write(6,810)imx(1),jmx(1),kmx(1), &
                imn(1),jmn(1),kmn(1)
           write(6,814)name,ndif,npos,xmx(1),xmn(1), &
                xmx(2),xmn(2)
           write(6,812)imx(2),jmx(2),kmx(2), &
                imn(2),jmn(2),kmn(2)

        else

           write(6,811)imx(1),jmx(1), &
                imn(1),jmn(1)
           write(6,814)name,ndif,npos,xmx(1),xmn(1), &
                xmx(2),xmn(2)
           write(6,813)imx(2),jmx(2), &
                imn(2),jmn(2)
        end if

        write(6,815) xbar(1), &
             xbar(2)

     end if              ! ndif.gt.0

     diffmw = sqrt(diffmw/wsum)
     write(6,'(a,a8,1pe11.4)') ' AREA WEIGHTED RMS ',name,diffmw

  else    ! Only 1 case

     if (numlev.gt.1) then                        ! 3-d field 
        write(6,810)imx(1),jmx(1),kmx(1), &
             imn(1),jmn(1),kmn(1)
        write(6,804)name,xmx(1),xmn(1)

     else                                         ! 2-d field

        write(6,811)imx(1),jmx(1), &
             imn(1),jmn(1)
        write(6,804)name,xmx(1),xmn(1)
     end if
     xbar(1) = xbar(1)/(nlon*nlat*numlev)
     write(6,805) xbar(1)
  end if

  return

803 format(1x,a8,1x,i6,1x,i6,1p2e23.15,e8.1,e23.15,e8.1,e23.15, &
       /,23x,               2e23.15,8x  ,e23.15,8x,  e23.15)
804 format(1x,a8,14x        ,1p2e23.15,e8.1,e23.15,e8.1,e23.15)
810 format(/24x,'(',i3,',',i3,',',i3,')', t47,'(',i3,',',i3,',',i3,')', &
       t79,'(',i3,',',i3,',',i3,')',t109,'(',i3,',',i3,',',i3,')')
812 format(24x ,'(',i3,',',i3,',',i3,')', t47,'(',i3,',',i3,',',i3,')')
811 format(/24x,'(',i3,',',i3,')'       , t47,'(',i3,',',i3,')', &
       t79,'(',i3,',',i3,')'       ,t109,'(',i3,',',i3,')')
813 format(24x ,'(',i3,',',i3,')'       , t47,'(',i3,',',i3,')')
805 format(10x,'avg abs field values:  ',1pe23.15,3x,'rms diff:',e8.1, &
       3x,'avg rel diff(npos): ',e8.1, &
       /,10x,'                       ',  e23.15,23x, &
       'avg decimal digits(ndif): ',0p,f4.1,' worst: ',f4.1)
814 format(1x,a8,1x,i6,1x,i6,1p2e23.15, &
       /23x,               2e23.15)
815 format(10x,'avg abs field values:  ',1p,e23.15, &
       /,10x,'                       ',   e23.15)

end subroutine printstats
