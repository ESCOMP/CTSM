module initHillslopeMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initializes (column-level) hillslope connectivity.
  !
  ! !USES:

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc,iam
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                

  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public initHillslopes  ! initialize hillslope connectivity
  public HillslopeDomPft ! change patch weights to a single pft
  !

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initHillslopes()
    !
    ! !DESCRIPTION: 
    ! Initialize hillslope connectivity.  Each landunit may have multiple
    ! hillslopes, and each hillslope may have multiple columns.  Specify
    ! uphill and downhill neighbors sequentially based on col%coli, col%colf.  
    ! Hilltop columns have no uphill neighbor; hillbottom columns have no 
    ! downhill neighbor.  

    !
    ! !USES
    use decompMod         , only : get_proc_bounds, get_clump_bounds, get_proc_clumps
    use clm_varctl        , only : nhillslope
    use landunit_varcon   , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: n,nc,nh,l,c,ci,cf    ! indices
    integer :: nclumps              ! number of clumps on this processor
    integer :: ncol_per_hill
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump
    integer :: begg, endg, begl, endl, begc, endc, begp, endp, &
         begCohort, endCohort
    !------------------------------------------------------------------------

    nclumps = get_proc_clumps()

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, l, nh, n, c)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)

       do l = bounds_clump%begl, bounds_clump%endl

! Set number of hillslopes; total number of columns already set 
! in clm_ptrs_compdown; this must be consistent with how columns set in 
! subgrid_get_info_* routines in subgridMod
! currently only _natveg specifies multiple columns
          if (lun%itype(l) == istsoil) then
             lun%nhillslopes(l) = nhillslope
          else
             lun%nhillslopes(l) = 0
          endif

          if (lun%nhillslopes(l) > 0) then

! number of columns per hillslope = lun%ncolumns / lun%nhillslopes 
! 1st column index of each hillslope is coli+(n-1)*ncol_per_hill for n(1:nhillslope)
! last column index is coli+n*ncol_per_hill-1
! This will be overwritten in HillslopeHydrologySurfaceDataMod
             ncol_per_hill = lun%ncolumns(l) / lun%nhillslopes(l)
             ! loop over hillslopes
             do nh = 1,lun%nhillslopes(l)
                ci = lun%coli(l) + (nh-1)*ncol_per_hill
                cf = lun%coli(l) + (nh)*ncol_per_hill - 1
                ! loop over columns within each hillslope
                do n = 1, ncol_per_hill
                   c = ci + (n-1)

                   col%hillslope_ndx(c) = nh

                   ! downhill columns (hillbottom has no downhill neighbor)
                   if(c < cf) then
                      col%cold(c) = c + 1
                   endif
                   ! uphill columns (hilltop has no uphill neighbor)
                   if(c > ci) then
                      col%colu(c) = c - 1
                   endif
                enddo ! end loop n
             enddo    ! end loop nh
          endif

       enddo ! end loop l
    enddo    ! end loop nc
    !$OMP END PARALLEL DO

  end subroutine initHillslopes

  !------------------------------------------------------------------------
  subroutine HillslopeDomPft()
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights such that each column has a single, 
    ! dominant pft.  

    !
    ! !USES
    use decompMod         , only : get_clump_bounds, get_proc_clumps
    use clm_varcon        , only : ispval
    use landunit_varcon   , only : istsoil, istcrop
    use PatchType         , only : patch
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: n,nc,p,pu,pl,l,c    ! indices
    integer :: nclumps             ! number of clumps on this processor
    integer :: upland_ivt  = 13 ! c3 non-arctic grass
    integer :: lowland_ivt = 8  ! broadleaf deciduous tree 
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump

    !------------------------------------------------------------------------

    nclumps = get_proc_clumps()

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, l, nh, n, c)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)

       do l = bounds_clump%begl, bounds_clump%endl

          if (lun%itype(l) == istsoil) then
             do c = lun%coli(l), lun%colf(l)

                sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
                sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
                sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))

                do p = col%patchi(c), col%patchf(c)
                   if(patch%itype(p) == lowland_ivt) pl = p
                   if(patch%itype(p) == upland_ivt)  pu = p
                enddo

                patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
                patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
                patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

                ! hillbottom
                if(col%cold(c) == ispval) then
                   patch%wtcol(pl) = sum_wtcol
                   patch%wtlunit(pl) = sum_wtlun
                   patch%wtgcell(pl) = sum_wtgrc
                else
                   patch%wtcol(pu) = sum_wtcol
                   patch%wtlunit(pu) = sum_wtlun
                   patch%wtgcell(pu) = sum_wtgrc
                endif
             enddo    ! end loop c
          endif
       enddo ! end loop l
    enddo    ! end loop nc
    !$OMP END PARALLEL DO

  end subroutine HillslopeDomPft

end module initHillslopeMod
