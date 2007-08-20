#!/bin/csh

set filename = spmdgs_subs.inc

touch $filename
mv -f $filename $filename.old
touch $filename

foreach DIM (1d 2d)
foreach TYPE (int real)

#  set filename = x${DIM}${TYPE}
#  rm -f $filename
#  touch $filename

  if ($DIM == 1d) then
     set CDIMS = "(:)    "
     set NDIMS = 1
     set VLB2  = 1
     set VUB2  = 1
     set V2NDL = ""
     set V2NDG = ""
  endif
  if ($DIM == 2d) then
     set CDIMS = "(:,:)  "
     set NDIMS = 2
     set VLB2  = "lbound(alocal,dim=2)"
     set VUB2  = "ubound(alocal,dim=2)"
     set V2NDL = ",n2"
     set V2NDG = ",n2-lb2+1"
  endif

  if ($TYPE == real) then
     set ATYPE  = "real(r8)"
     set NCTYPE = "double"
     set MPTYPE = "MPI_REAL8"
     set MCLIST = "rList"
     set MCATTR = "Rattr"
     set MCSTR  = "rstring"
  endif
  if ($TYPE == int) then
     set ATYPE  = "integer "
     set NCTYPE = "int"
     set MPTYPE = "MPI_INTEGER"
     set MCLIST = "iList"
     set MCATTR = "Iattr"
     set MCSTR  = "istring"
  endif

echo $DIM $TYPE

cat >> $filename <<EOF

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_${DIM}array_${TYPE}
!
! !INTERFACE:
  subroutine scatter_${DIM}array_${TYPE} (alocal, aglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter ${TYPE} ${DIM} array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    ${ATYPE}, pointer            :: alocal${CDIMS}   ! local  data (output)
    ${ATYPE}, pointer            :: aglobal${CDIMS}  ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    ${ATYPE},pointer   :: adata(:)   ! local data array
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    integer, pointer,dimension(:) :: perm    ! mct permuter
    character(len=*),parameter :: subname = 'scatter_${DIM}array_${TYPE}'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')
    call get_clmlevel_gsmap(clmlevel,gsmap,perm)

    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = ${VLB2}
    ub2 = ${VUB2}

    rstring = ""
    istring = ""

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(${MCSTR}) == 0) then
          ${MCSTR} = trim(fname)
       else
          ${MCSTR} = trim(${MCSTR})//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    if (debug > 1) call t_startf(trim(subname)//'_pack')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          adata(1:lsize) = aglobal(1:lsize${V2NDG})
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_import${MCATTR}(AVi,trim(fname),adata,lsize)
       enddo
       deallocate(adata)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_scat')

    call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)
    call mct_aVect_unpermute(AVo, perm, dieWith=subname)

    if (debug > 1) call t_stopf(trim(subname)//'_scat')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    lsize = size(alocal,dim=1)
    allocate(adata(lsize))
    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_export${MCATTR}(AVo,trim(fname),adata,lsize)
       do n1 = lb1,ub1
          alocal(n1${V2NDL}) = adata(n1-lb1+1)
       enddo
    enddo
    deallocate(adata)

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVi)
    endif
    call mct_aVect_clean(AVo)

    call t_stopf(trim(subname)//'_total')

  end subroutine scatter_${DIM}array_${TYPE}

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_${DIM}array_${TYPE}
!
! !INTERFACE:
  subroutine gather_${DIM}array_${TYPE} (alocal, aglobal, clmlevel, missing)
!
! !DESCRIPTION:
! Wrapper routine to gather ${TYPE} ${DIM} array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    ${ATYPE}, pointer            :: alocal${CDIMS}   ! local  data (output)
    ${ATYPE}, pointer            :: aglobal${CDIMS}  ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
    ${ATYPE},optional,intent(in) :: missing     ! missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    ${ATYPE},pointer   :: adata(:)   ! temporary data array
    integer ,pointer   :: mvect(:)   ! local array for mask
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    integer, pointer,dimension(:) :: perm    ! mct permuter
    character(len=*),parameter :: subname = 'gather_${DIM}array_${TYPE}'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')

    call get_clmlevel_gsmap(clmlevel,gsmap,perm)

    lsize = size(alocal,dim=1)
    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = ${VLB2}
    ub2 = ${VUB2}
   
    rstring = ""
    istring = ""

    if (present(missing)) then
       istring = "mask"
    endif

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(${MCSTR}) == 0) then
          ${MCSTR} = trim(fname)
       else
          ${MCSTR} = trim(${MCSTR})//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)

    if (debug > 1) call t_startf(trim(subname)//'_pack')
    allocate(adata(lsize))
    do n2 = lb2,ub2
       do n1 = lb1,ub1
          adata(n1-lb1+1) = alocal(n1${V2NDL})
       enddo
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_import${MCATTR}(AVi,trim(fname),adata,lsize)
    enddo
    deallocate(adata)

    if (present(missing)) then
       allocate(mvect(lsize))
       do n1 = lb1,ub1
          mvect(n1-lb1+1) = 1
       enddo
       call mct_aVect_importIattr(AVi,"mask",mvect,lsize)
       deallocate(mvect)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_gath')

    call mct_aVect_permute(AVi, perm, dieWith=subname)
    if (present(missing)) then
! tcx wait for update in mct, then get rid of "mask"
!       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom, missing = missing)
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    else
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_gath')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_export${MCATTR}(AVo,trim(fname),adata,lsize)
          aglobal(1:lsize${V2NDG}) = adata(1:lsize)
       enddo
       deallocate(adata)
       if (present(missing)) then
          allocate(mvect(lsize))
          call mct_aVect_exportIattr(AVo,"mask",mvect,lsize)
          do n1 = 1,lsize
             if (mvect(n1) == 0) then
                do n2 = lb2,ub2
                   aglobal(n1${V2NDG}) = missing
                enddo
             endif
          enddo
          deallocate(mvect)
       endif
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVo)
    endif

    call mct_aVect_clean(AVi)

    call t_stopf(trim(subname)//'_total')

  end subroutine gather_${DIM}array_${TYPE}

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allgather_${DIM}array_${TYPE}
!
! !INTERFACE:
  subroutine allgather_${DIM}array_${TYPE} (alocal, aglobal, clmlevel, missing)
!
! !DESCRIPTION:
! Wrapper routine to perform an allgatherv of ${DIM} ${TYPE} array
!
! !ARGUMENTS:
    implicit none
    ${ATYPE}, pointer            :: alocal${CDIMS}   ! local  data (output)
    ${ATYPE}, pointer            :: aglobal${CDIMS}  ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
    ${ATYPE},optional,intent(in) :: missing     ! missing value
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier          ! error code
    character(len=*),parameter :: subname = 'allgather_${DIM}array_${TYPE}'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname)
    endif

    if (present(missing)) then
       call gather_data_to_master(alocal,aglobal,clmlevel,missing)
    else
       call gather_data_to_master(alocal,aglobal,clmlevel)
    endif
    call mpi_bcast (aglobal, size(aglobal), ${MPTYPE}, 0, mpicom, ier)
    if (ier/=0 ) then
       write(iulog,*) trim(subname),ier
       call endrun()
    endif
    call t_stopf(trim(subname)//'_total')

  end subroutine allgather_${DIM}array_${TYPE}


EOF

end
end

exit
