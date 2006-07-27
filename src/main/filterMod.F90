#include <misc.h>
#include <preproc.h>

module filterMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: filterMod
!
! !DESCRIPTION:
! Module of filters used for processing columns and pfts of particular
! types, including lake, non-lake, soil, snow, non-snow, and
! naturally-vegetated patches.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save

  type clumpfilter
#ifdef DGVM
     integer, pointer :: natvegp(:) ! DGVM naturally-vegetated (present)
                                    ! filter (pfts)
     integer :: num_natvegp         ! number of pfts in naturally-vegetated
                                    ! filter
#endif
     integer, pointer :: lakep(:)   ! lake filter (pfts)
     integer :: num_lakep           ! number of pfts in lake filter
     integer, pointer :: nolakep(:) ! non-lake filter (pfts)
     integer :: num_nolakep         ! number of pfts in non-lake filter
     integer, pointer :: lakec(:)   ! lake filter (columns)
     integer :: num_lakec           ! number of columns in lake filter
     integer, pointer :: nolakec(:) ! non-lake filter (columns)
     integer :: num_nolakec         ! number of columns in non-lake filter
     integer, pointer :: soilc(:)   ! soil filter (columns)
     integer :: num_soilc           ! number of columns in soil filter 
     integer, pointer :: soilp(:)   ! soil filter (pfts)
     integer :: num_soilp           ! number of pfts in soil filter 
     integer, pointer :: snowc(:)   ! snow filter (columns) 
     integer :: num_snowc           ! number of columns in snow filter 
     integer, pointer :: nosnowc(:) ! non-snow filter (columns) 
     integer :: num_nosnowc         ! number of columns in non-snow filter 
  end type clumpfilter
  type(clumpfilter), allocatable, public :: filter(:)
!
  public allocFilters   ! allocate memory for filters
  public setFilters     ! set filters
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 11/13/03, Peter Thornton: Added soilp and num_soilp
!
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocFilters
!
! !INTERFACE:
  subroutine allocFilters()
!
! !DESCRIPTION:
! Allocate CLM filters.
!
! !USES:
    use clmtype
    use decompMod , only : get_proc_clumps, get_clump_bounds
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2004.04.27 DGVM naturally-vegetated filter added by Forrest Hoffman
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: nc          ! clump index
    integer :: nclumps     ! total number of clumps on this processor
    integer :: begp, endp  ! per-clump beginning and ending pft indices
    integer :: begc, endc  ! per-clump beginning and ending column indices
    integer :: begl, endl  ! per-clump beginning and ending landunit indices
    integer :: begg, endg  ! per-clump beginning and ending gridcell indices
    integer :: ier         ! error status
!------------------------------------------------------------------------

    ! Determine clump variables for this processor

    nclumps = get_proc_clumps()
    ier = 0
    if( .not. allocated(filter)) then
       allocate(filter(nclumps), stat=ier)
    end if
    if (ier /= 0) then
       write (6,*) 'allocFilters(): allocation error for clumpsfilters'
       call endrun
    end if

    ! Loop over clumps on this processor

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
!WARNING: Do not put CSDs around loops which call allocate()
    do nc = 1, nclumps
       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
       allocate(filter(nc)%lakec(endc-begc+1))
       allocate(filter(nc)%nolakec(endc-begc+1))
       allocate(filter(nc)%lakep(endp-begp+1))
       allocate(filter(nc)%nolakep(endp-begp+1))
       allocate(filter(nc)%soilc(endc-begc+1))
       allocate(filter(nc)%soilp(endp-begp+1))
       allocate(filter(nc)%snowc(endc-begc+1))
       allocate(filter(nc)%nosnowc(endc-begc+1))
#ifdef DGVM
       allocate(filter(nc)%natvegp(endp-begp+1))
#endif
    end do
!$OMP END PARALLEL DO

  end subroutine allocFilters

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setFilters
!
! !INTERFACE:
  subroutine setFilters()
!
! !DESCRIPTION:
! Set CLM filters.
!
! !USES:
    use clmtype
    use decompMod , only : get_proc_clumps, get_clump_bounds
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2004.04.27 DGVM naturally-vegetated filter added by Forrest Hoffman
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: nc          ! clump index
    integer :: c,l,p       ! column, landunit, pft indices
    integer :: nclumps     ! total number of clumps on this processor
    integer :: fl          ! lake filter index
    integer :: fnl         ! non-lake filter index
    integer :: fs          ! soil filter index
    integer :: begp, endp  ! per-clump beginning and ending pft indices
    integer :: begc, endc  ! per-clump beginning and ending column indices
    integer :: begl, endl  ! per-clump beginning and ending landunit indices
    integer :: begg, endg  ! per-clump beginning and ending gridcell indices
!------------------------------------------------------------------------

    ! Loop over clumps on this processor

    nclumps = get_proc_clumps()
!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp,fl,fnl,fs,p,c,l)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp,fl,fnl,fs,p,c,l)
#endif
    do nc = 1,nclumps

       ! Determine clump boundaries

       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

       ! Create lake and non-lake filters at column-level 

       fl = 0
       fnl = 0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = clm3%g%l%c%landunit(c)
          if (clm3%g%l%lakpoi(l)) then
             fl = fl + 1
             filter(nc)%lakec(fl) = c
          else
             fnl = fnl + 1
             filter(nc)%nolakec(fnl) = c
          end if
       end do
       filter(nc)%num_lakec = fl
       filter(nc)%num_nolakec = fnl

       ! Create lake and non-lake filters at pft-level 
       ! Filter will only be active if weight of pft wrt gcell is nonzero

       fl = 0
       fnl = 0
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (clm3%g%l%c%p%wtgcell(p) > 0._r8) then
             l = clm3%g%l%c%p%landunit(p)
             if (clm3%g%l%lakpoi(l) ) then
                fl = fl + 1
                filter(nc)%lakep(fl) = p
             else
                fnl = fnl + 1
                filter(nc)%nolakep(fnl) = p
             end if
          end if
       end do
       filter(nc)%num_lakep = fl
       filter(nc)%num_nolakep = fnl

       ! Create soil filter at column-level

       fs = 0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = clm3%g%l%c%landunit(c)
          if (clm3%g%l%itype(l) == istsoil) then
             fs = fs + 1
             filter(nc)%soilc(fs) = c
          end if
       end do
       filter(nc)%num_soilc = fs

       ! Create soil filter at pft-level
       ! Filter will only be active if weight of pft wrt gcell is nonzero

       fs = 0
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (clm3%g%l%c%p%wtgcell(p) > 0._r8) then
             l = clm3%g%l%c%p%landunit(p)
             if (clm3%g%l%itype(l) == istsoil) then
                fs = fs + 1
                filter(nc)%soilp(fs) = p
             end if
          end if
       end do
       filter(nc)%num_soilp = fs

       ! Note: snow filters are reconstructed each time step in Hydrology2
       ! Note: DGVM present vegetated filter is reconstructed each time DGVM is run

    end do
!$OMP END PARALLEL DO
#if !defined (USE_OMP)
!CSD$ END PARALLEL DO
#endif

  end subroutine setFilters

end module filterMod
