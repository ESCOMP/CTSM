module clm_glclnd

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_glclnd
!
! !DESCRIPTION:
! Handle arrays used for exchanging data between glc and land model.
! Based on clm_atmlnd (but without mapping routines because glc data
!  is send and received on the lnd decomposition, at least for now).
!
! The fields sent from the lnd component to the glc component via
!  the coupler are labeled 's2x', or sno to coupler.
! The fields received by the lnd component from the glc component
!  via the coupler are labeled 'x2s', or coupler to sno.
! 'Sno' is a misnomer in that the exchanged data are related to
!  the ice beneath the snow, not the snow itself.  But by CESM convention,
! 'ice' refers to sea ice, not land ice.
!
! !USES:
  use decompMod   , only : get_proc_bounds, get_proc_bounds_atm
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use spmdMod     , only : masterproc
  use clm_varctl  , only : iulog
  use clm_varctl  , only : glc_nec
!
! !REVISION HISTORY:
! Created by William Lipscomb, Dec. 2007, based on clm_atmlnd.F90.
!
! !PUBLIC TYPES:
  implicit none

!----------------------------------------------------
! glc -> land variables structure
!----------------------------------------------------
  type glc2lnd_type
     real(r8), pointer :: frac(:,:) 
     real(r8), pointer :: topo(:,:) 
     real(r8), pointer :: rofi(:,:) 
     real(r8), pointer :: rofl(:,:) 
     real(r8), pointer :: hflx(:,:) 
  end type glc2lnd_type

!----------------------------------------------------
! land -> glc variables structure
!----------------------------------------------------
  type lnd2glc_type
     real(r8), pointer :: tsrf(:,:) 
     real(r8), pointer :: topo(:,:)
     real(r8), pointer :: qice(:,:)
  end type lnd2glc_type

  type (lnd2glc_type), public, target :: clm_s2x  ! s2x fields on clm grid
  type (glc2lnd_type), public, target :: clm_x2s  ! x2s fields on clm grid

  type (lnd2glc_type), public, target :: atm_s2x  ! s2x fields on atm grid
  type (glc2lnd_type), public, target :: atm_x2s  ! x2s fields on atm grid

! !PUBLIC MEMBER FUNCTIONS:
  public :: init_glc2lnd_type
  public :: init_lnd2glc_type
  public :: clm_maps2x
  public :: clm_mapx2s
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!----------------------------------------------------

contains


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_glc2lnd_type
!
! !INTERFACE:
  subroutine init_glc2lnd_type(beg, end, x2s)
!
! !DESCRIPTION:
! Initialize glc variables required by the land
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (glc2lnd_type), intent(inout):: x2s
!
! !REVISION HISTORY:
! Created by William Lipscomb, based on init_atm2lnd_type
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(x2s%frac(beg:end,glc_nec))
  allocate(x2s%topo(beg:end,glc_nec))
  allocate(x2s%rofi(beg:end,glc_nec))
  allocate(x2s%rofl(beg:end,glc_nec))
  allocate(x2s%hflx(beg:end,glc_nec))

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  x2s%frac(beg:end,:) = ival
  x2s%topo(beg:end,:) = ival
  x2s%rofi(beg:end,:) = ival
  x2s%rofl(beg:end,:) = ival
  x2s%hflx(beg:end,:) = ival

end subroutine init_glc2lnd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2glc_type
!
! !INTERFACE:
  subroutine init_lnd2glc_type(beg, end, s2x)
!
! !DESCRIPTION:
! Initialize land variables required by glc
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (lnd2glc_type), intent(inout):: s2x
!
! !REVISION HISTORY:
! Created by William Lipscomb, based on init_lnd2atm_type
!
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(s2x%tsrf(beg:end,glc_nec))
  allocate(s2x%topo(beg:end,glc_nec))
  allocate(s2x%qice(beg:end,glc_nec))

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  s2x%tsrf(beg:end,:) = ival
  s2x%topo(beg:end,:) = ival
  s2x%qice(beg:end,:) = ival

end subroutine init_lnd2glc_type

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_maps2x
!
! !INTERFACE:
  subroutine clm_maps2x(s2x_src, s2x_dst)
!
! !DESCRIPTION:
! Maps lnd2glc fields from clm grid to external grid
!
! !USES:
  use downscaleMod  , only : map_maparrayl, map1dl_l2a
!
! !ARGUMENTS:
  implicit none
  type(lnd2glc_type), intent(in)  :: s2x_src
  type(lnd2glc_type), intent(out) :: s2x_dst
!
! !REVISION HISTORY:
! Created by William Lipscomb based on clm_mapl2a 
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
!------------------------------------------------------------------------------
 
  !--- allocate temporaries
  call get_proc_bounds    (begg_s, endg_s)
  call get_proc_bounds_atm(begg_d, endg_d)
 
  nflds = 3    
 
  allocate(asrc(begg_s:endg_s,nflds))
  allocate(adst(begg_d:endg_d,nflds))

  do n = 1, glc_nec
 
     ix = 0
     ix=ix+1; asrc(:,ix) = s2x_src%tsrf(:,n)  
     ix=ix+1; asrc(:,ix) = s2x_src%topo(:,n)  
     ix=ix+1; asrc(:,ix) = s2x_src%qice(:,n)  
 
     call map_maparrayl(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, map1dl_l2a)
 
     ix = 0
     ix=ix+1; s2x_dst%tsrf(:,n) = adst(:,ix)
     ix=ix+1; s2x_dst%topo(:,n) = adst(:,ix)
     ix=ix+1; s2x_dst%qice(:,n) = adst(:,ix)
 
  enddo  ! loop over elevation classes

  deallocate(asrc)
  deallocate(adst)
 
end subroutine clm_maps2x
 
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapx2s
!
! !INTERFACE:
  subroutine clm_mapx2s(x2s_src, x2s_dst)
!
! !DESCRIPTION:
! Maps glc2lnd fields from external grid to clm grid
!
! !USES:
  use downscaleMod, only : map_maparrayl, map1dl_a2l, map1dl_l2a, map_setptrs
  use decompMod   , only : ldecomp,adecomp
  use domainMod   , only : ldomain,adomain
!
! !ARGUMENTS:
  implicit none
  type(glc2lnd_type), intent(in)  :: x2s_src
  type(glc2lnd_type), intent(out) :: x2s_dst
!
! !REVISION HISTORY:
! Created by William Lipscomb based on clm_mapa2l 
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
  logical  :: first_call = .true.
!------------------------------------------------------------------------------
 
  if (first_call .and. masterproc) then
    write(iulog,*) 'clm_mapx2s subroutine'
  endif
 
  !--- allocate temporaries
  call get_proc_bounds_atm(begg_s, endg_s)
  call get_proc_bounds    (begg_d, endg_d)
 
  nflds = 5
 
  allocate(asrc(begg_s:endg_s,nflds))
  allocate(adst(begg_d:endg_d,nflds))

  do n = 1, glc_nec
 
     ix = 0
     ix=ix+1; asrc(:,ix) = x2s_src%frac(:,n)  
     ix=ix+1; asrc(:,ix) = x2s_src%topo(:,n)  
     ix=ix+1; asrc(:,ix) = x2s_src%rofi(:,n)  
     ix=ix+1; asrc(:,ix) = x2s_src%rofl(:,n)  
     ix=ix+1; asrc(:,ix) = x2s_src%hflx(:,n)  
 
     call map_maparrayl(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, map1dl_a2l)
 
     ix = 0
     ix=ix+1; x2s_dst%frac(:,n)  =   adst(:,ix)
     ix=ix+1; x2s_dst%topo(:,n)  =   adst(:,ix)
     ix=ix+1; x2s_dst%rofi(:,n)  =   adst(:,ix)
     ix=ix+1; x2s_dst%rofl(:,n)  =   adst(:,ix)
     ix=ix+1; x2s_dst%hflx(:,n)  =   adst(:,ix)

  enddo

  deallocate(asrc)
  deallocate(adst)

  if (first_call.and.masterproc) then
    write(iulog,*) 'clm_mapx2s mapping complete'
  endif

  first_call = .false.

  end subroutine clm_mapx2s
 
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_clm_s2x
!
! !INTERFACE:
  subroutine create_clm_s2x(clm_s2x)
!
! !DESCRIPTION:
! Assign values to clm_s2x based on the appropriate derived types
!
! !USES:
  use clm_varctl  , only : glc_smb
  use clmtype     , only : clm3
  use domainMod   , only : ldomain
  use clm_varcon  , only : istice_mec
  use clm_atmlnd  , only : clm_l2a, clm_a2l
  use clm_varcon  , only : spval
!
! !ARGUMENTS:
  implicit none

  type(lnd2glc_type), intent(out) :: clm_s2x
!
! !REVISION HISTORY:
! Written by William Lipscomb, Feb. 2009 
!

    integer :: begg, endg              ! per-proc beginning and ending gridcell indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: c, l, g, n              ! indices
    integer , pointer :: ityplun(:)    ! landunit type
    integer , pointer :: clandunit(:)  ! column's landunit index
    integer , pointer :: cgridcell(:)  ! column's gridcell index

    ! Assign local pointers to derived type members

    clandunit     => clm3%g%l%c%landunit
    cgridcell     => clm3%g%l%c%gridcell
    ityplun       => clm3%g%l%itype

    ! Get processor bounds

    call get_proc_bounds(begg, endg, begc=begc, endc=endc)

    ! initialize to be safe
    clm_s2x%tsrf(:,:) = 0._r8
    clm_s2x%topo(:,:) = 0._r8
    clm_s2x%qice(:,:) = 0._r8

    ! fill the clm_s2x vector on the clm grid

    if (glc_smb) then   ! send surface mass balance info

       do c = begc, endc
          l = clandunit(c)
          g = cgridcell(c)

          if (ityplun(l) == istice_mec) then
             n = c - clm3%g%l%coli(l) + 1    ! elevation class index
                                             ! (assumes all elevation classes are populated)
             clm_s2x%tsrf(g,n) = clm3%g%l%c%ces%t_soisno(c,1)
             clm_s2x%qice(g,n) = clm3%g%l%c%cwf%qflx_glcice(c)
             clm_s2x%topo(g,n) = clm3%g%l%c%cps%glc_topo(c)

             ! Check for bad values of qice
             if ( abs(clm_s2x%qice(g,n)) > 1.0_r8 .and. clm_s2x%qice(g,n) /= spval) then
                write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, clm_s2x%qice(g,n)
             endif

           endif    ! istice_mec
       enddo        ! c

    else  ! Pass PDD info (same info in each elevation class)
             ! It might make sense to require glc_nec = 1 for this case
             

       do n = 1, glc_nec
          do g = begg, endg
             clm_s2x%tsrf(g,n) = clm_l2a%t_ref2m(g)
             clm_s2x%qice(g,n) = clm_a2l%forc_snow(g)   ! Assume rain runs off
             clm_s2x%topo(g,n) = ldomain%topo(g)

             ! Check for bad values of qice
             if (clm_s2x%qice(g,n) > -1.0_r8 .and. clm_s2x%qice(g,n) < 1.0_r8) then
                continue
             else
                write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, clm_s2x%qice(g,n)
                write(iulog,*) 'forc_rain =', clm_a2l%forc_rain(g)
                write(iulog,*) 'forc_snow =', clm_a2l%forc_snow(g)
             endif

          enddo
       enddo

    endif   ! glc_smb

end subroutine create_clm_s2x

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: unpack_clm_x2s
!
! !INTERFACE:
  subroutine unpack_clm_x2s(clm_x2s)
!
! !DESCRIPTION:
! Unpack clm_x2s and update the appropriate derived types
!
! !USES:
  use clm_varcon  , only : istice_mec
  use clmtype     , only : clm3

!
! !ARGUMENTS:
  implicit none

  type(glc2lnd_type), intent(in) :: clm_x2s
!
! !REVISION HISTORY:
! Written by William Lipscomb, Feb. 2009 
!
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: c, l, g, n              ! indices
    integer , pointer :: ityplun(:)    ! landunit type
    integer , pointer :: clandunit(:)  ! column's landunit index
    integer , pointer :: cgridcell(:)  ! column's gridcell index

    logical :: update_glc2sno_fields   ! if true, update glacier_mec fields

    ! Assign local pointers to derived type members

    clandunit     => clm3%g%l%c%landunit
    cgridcell     => clm3%g%l%c%gridcell
    ityplun       => clm3%g%l%itype

    update_glc2sno_fields = .false.
                                           
    if (update_glc2sno_fields) then 
   
       do c = begc, endc
          l = clandunit(c)
          g = cgridcell(c)

          if (ityplun(l) == istice_mec) then
             n = c - clm3%g%l%coli(l) + 1    ! elevation class index
             clm3%g%l%c%cps%glc_frac(c) = clm_x2s%frac(g,n)
             clm3%g%l%c%cps%glc_topo(c) = clm_x2s%topo(g,n)
             clm3%g%l%c%cwf%glc_rofi(c) = clm_x2s%rofi(g,n)
             clm3%g%l%c%cwf%glc_rofl(c) = clm_x2s%rofl(g,n)
             clm3%g%l%c%cef%eflx_bot(c) = clm_x2s%hflx(g,n)

          endif
       enddo

    endif   ! update fields

    end subroutine unpack_clm_x2s

!------------------------------------------------------------------------

end module clm_glclnd

