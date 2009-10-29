#include <misc.h>
#include <preproc.h>

module iobinary

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: iobinary
!
! !DESCRIPTION:
! Set of wrappers to write binary I/O
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod        , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER, &
                              MPI_LOGICAL
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
  use decompMod      , only : get_clmlevel_gsize
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  interface bin_iolocal
     module procedure bin_2darray_int
     module procedure bin_2darray_real
  end interface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Updated by tcraig, 3/2007
!
!
! !PRIVATE MEMBER FUNCTIONS: None
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: bin_2darray_int
!
! !INTERFACE:
  subroutine bin_2darray_int (iu, arrayin, clmlevel, flag)
!
! !DESCRIPTION:
! Wrapper routine to read/write integer 2d array from restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    integer, pointer    :: arrayin(:,:)        !input data
    character(len=*), intent(in) :: clmlevel   !type of input data
    character(len=*), intent(in) :: flag       !'read' or 'write'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Updated by tcraig, 3/2007
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: n                               !index
    integer :: ier                             !return code
    integer :: gsize                           !size of first dimension
    integer :: lb1,ub1,lb2,ub2                 !bound of arrayin
    integer, pointer :: arrayl(:)              !temporary
    integer, pointer :: arrayg(:)              !temporary
    character(len=*),parameter :: subname = 'bin_2darray_int'
!-----------------------------------------------------------------------

    if (flag /= 'read' .and. flag /= 'write') then
       write(iulog,*) trim(subname),' error in flag ',trim(flag)
       call endrun()
    endif

    gsize = get_clmlevel_gsize(clmlevel)
    if (masterproc) then
       allocate(arrayg(gsize),stat=ier)
       if (ier /= 0) then
          write(iulog,*) trim(subname),'arrayg allocation error'
          call endrun()
       end if
    endif

    lb1 = lbound(arrayin,dim=1)
    ub1 = ubound(arrayin,dim=1)
    lb2 = lbound(arrayin,dim=2)
    ub2 = ubound(arrayin,dim=2)

    allocate(arrayl(lb1:ub1),stat=ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname),'arrayg allocation error'
       call endrun()
    end if

    do n = lb2,ub2

       if (flag == 'write') then
          arrayl(lb1:ub1) = arrayin(lb1:ub1,n)
          call gather_data_to_master(arrayl, arrayg, clmlevel)
       endif

      if (masterproc) then
         if (flag == 'write') write (iu,iostat=ier) arrayg
         if (flag == 'read' ) read  (iu,iostat=ier) arrayg
         if (ier /= 0 ) then
            write(iulog,*) trim(subname),'ier = ',ier,' on i/o unit = ',iu
            call endrun()
         endif
      endif

      if (flag == 'read') then
          call scatter_data_from_master(arrayl, arrayg, clmlevel)
          arrayin(lb1:ub1,n) = arrayl(lb1:ub1)
      endif

    enddo

    if (masterproc) deallocate(arrayg)
    deallocate(arrayl)

  end subroutine bin_2darray_int


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: bin_2darray_real
!
! !INTERFACE:
  subroutine bin_2darray_real (iu, arrayin, clmlevel, flag)
!
! !DESCRIPTION:
! Wrapper routine to read/write integer 2d array from restart binary file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iu                  !input unit
    real(r8), pointer   :: arrayin(:,:)        !input data
    character(len=*), intent(in) :: clmlevel   !type of input data
    character(len=*), intent(in) :: flag       !'read' or 'write'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Updated by tcraig, 3/2007
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: n                               !index
    integer :: ier                             !return code
    integer :: gsize                           !size of first dimension
    integer :: lb1,ub1,lb2,ub2                 !bound of arrayin
    real(r8), pointer :: arrayl(:)             !temporary
    real(r8), pointer :: arrayg(:)             !temporary
    character(len=*),parameter :: subname = 'bin_2darray_real'
!-----------------------------------------------------------------------

    if (flag /= 'read' .and. flag /= 'write') then
       write(iulog,*) trim(subname),' error in flag ',trim(flag)
       call endrun()
    endif

    gsize = get_clmlevel_gsize(clmlevel)
    if (masterproc) then
       allocate(arrayg(gsize),stat=ier)
       if (ier /= 0) then
          write(iulog,*) trim(subname),'arrayg allocation error'
          call endrun()
       end if
    endif

    lb1 = lbound(arrayin,dim=1)
    ub1 = ubound(arrayin,dim=1)
    lb2 = lbound(arrayin,dim=2)
    ub2 = ubound(arrayin,dim=2)

    allocate(arrayl(lb1:ub1),stat=ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname),'arrayg allocation error'
       call endrun()
    end if

    do n = lb2,ub2

       if (flag == 'write') then
          arrayl(lb1:ub1) = arrayin(lb1:ub1,n)
          call gather_data_to_master(arrayl, arrayg, clmlevel)
       endif

      if (masterproc) then
         if (flag == 'write') write (iu,iostat=ier) arrayg
         if (flag == 'read' ) read  (iu,iostat=ier) arrayg
         if (ier /= 0 ) then
            write(iulog,*) trim(subname),'ier = ',ier,' on i/o unit = ',iu
            call endrun()
         endif
      endif

      if (flag == 'read') then
          call scatter_data_from_master(arrayl, arrayg, clmlevel)
          arrayin(lb1:ub1,n) = arrayl(lb1:ub1)
      endif

    enddo

    if (masterproc) deallocate(arrayg)
    deallocate(arrayl)

  end subroutine bin_2darray_real


end module iobinary
