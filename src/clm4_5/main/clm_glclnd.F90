module clm_glclnd

  !-----------------------------------------------------------------------
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
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use decompMod   , only : get_proc_bounds
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use clm_varpar  , only : maxpatch_glcmec
  use clm_varctl  , only : iulog, glc_smb
  use abortutils  , only : endrun
  use decompMod   , only : bounds_type
  !
  ! !REVISION HISTORY:
  ! Created by William Lipscomb, Dec. 2007, based on clm_atmlnd.F90.
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! glc -> land variables structure
  type, public :: glc2lnd_type
     real(r8), pointer :: frac(:,:) => null()
     real(r8), pointer :: topo(:,:) => null()
     real(r8), pointer :: hflx(:,:) => null()
  end type glc2lnd_type

  ! land -> glc variables structure
  type, public :: lnd2glc_type
     real(r8), pointer :: tsrf(:,:) => null()
     real(r8), pointer :: topo(:,:) => null()
     real(r8), pointer :: qice(:,:) => null()
  end type lnd2glc_type

  type (lnd2glc_type), public, target :: clm_s2x  ! s2x fields on clm grid
  type (glc2lnd_type), public, target :: clm_x2s  ! x2s fields on clm grid

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_glc2lnd_type
  public :: init_lnd2glc_type
  public :: update_clm_s2x
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_glc2lnd_type(bounds, x2s)
    !
    ! !DESCRIPTION:
    ! Initialize glc variables required by the land
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    type (glc2lnd_type), intent(inout):: x2s
    !------------------------------------------------------------------------

    allocate(x2s%frac(bounds%begg:bounds%endg,maxpatch_glcmec))
    x2s%frac(bounds%begg:bounds%endg,maxpatch_glcmec)=0.0_r8
    allocate(x2s%topo(bounds%begg:bounds%endg,maxpatch_glcmec))
    x2s%topo(bounds%begg:bounds%endg,maxpatch_glcmec)=0.0_r8
    allocate(x2s%hflx(bounds%begg:bounds%endg,maxpatch_glcmec))
    x2s%hflx(bounds%begg:bounds%endg,maxpatch_glcmec)=0.0_r8

  end subroutine init_glc2lnd_type

  !------------------------------------------------------------------------
  subroutine init_lnd2glc_type(bounds, s2x)
    !
    ! !DESCRIPTION:
    ! Initialize land variables required by glc
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    type (lnd2glc_type), intent(inout):: s2x
    !------------------------------------------------------------------------

    allocate(s2x%tsrf(bounds%begg:bounds%endg,maxpatch_glcmec))
    s2x%tsrf(bounds%begg:bounds%endg,maxpatch_glcmec)=0.0_r8
    allocate(s2x%topo(bounds%begg:bounds%endg,maxpatch_glcmec))
    s2x%topo(bounds%begg:bounds%endg,maxpatch_glcmec)=0.0_r8
    allocate(s2x%qice(bounds%begg:bounds%endg,maxpatch_glcmec))
    s2x%qice(bounds%begg:bounds%endg,maxpatch_glcmec)=0.0_r8

  end subroutine init_lnd2glc_type

  !------------------------------------------------------------------------------
  subroutine update_clm_s2x(bounds, init)
    !
    ! !DESCRIPTION:
    ! Assign values to clm_s2x based on the appropriate derived types
    !
    ! !USES:
    use clmtype 
    use domainMod , only : ldomain
    use clm_varcon, only : istice_mec
    use clm_atmlnd, only : clm_l2a, clm_a2l
    use clm_varcon, only : spval
    !
    ! !ARGUMENTS:
    implicit none
    logical, intent(in)   :: init    ! if true=>only set a subset of fields
    type(bounds_type), optional, intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begg, endg              ! per-proc beginning and ending gridcell indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: c, l, g, n              ! indices
    !------------------------------------------------------------------------------

    if (.not. present(bounds)) then 
       call get_proc_bounds(begg=begg, endg=endg, begc=begc, endc=endc)
    end if

    ! Initialize qice because otherwise it will remain unset if init=true and
    ! glc_smb=true; note that the value here is the value qice will keep if these
    ! conditions hold

    clm_s2x%qice(:,:) = 0._r8

    ! and initialize the other variables just to be safe

    clm_s2x%tsrf(:,:) = 0._r8
    clm_s2x%topo(:,:) = 0._r8

    ! Fill the clm_s2x vector on the clm grid

    if (glc_smb) then   ! send surface mass balance info
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          g = col%gridcell(c)

          ! Following assumes all elevation classes are populated
          if (lun%itype(l) == istice_mec) then
             n = c -lun%coli(l) + 1    ! elevation class index

             ! t_soisno and glc_topo are valid even in initialization, so tsrf and topo
             ! are set here regardless of the value of init. But qflx_glcice is not valid
             ! until the run loop; thus, in initialization, we will use the default value
             ! for qice, as set above.
             clm_s2x%tsrf(g,n) = ces%t_soisno(c,1)
             clm_s2x%topo(g,n) = cps%glc_topo(c)
             if (.not. init) then
                clm_s2x%qice(g,n) = cwf%qflx_glcice(c)

                ! Check for bad values of qice
                if ( abs(clm_s2x%qice(g,n)) > 1.0_r8 .and. clm_s2x%qice(g,n) /= spval) then
                   write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, clm_s2x%qice(g,n)
                endif
             end if

           endif    ! istice_mec
       enddo        ! c
    else  ! Pass PDD info (same info in each elevation class)
       ! Require maxpatch_glcmec = 1 for this case 
       if (maxpatch_glcmec .ne. 1) then
          call endrun('update_clm_s2x error: maxpatch_glcmec must be 1 if glc_smb is false') 
       end if
       n = 1
       do g = bounds%begg,bounds%endg
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
    endif   ! glc_smb

  end subroutine update_clm_s2x

end module clm_glclnd

