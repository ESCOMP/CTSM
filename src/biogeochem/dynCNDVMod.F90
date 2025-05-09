module dynCNDVMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle weight updates associated with prognostic dynamic vegetation (CNDV)
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use decompMod    , only : bounds_type
  use LandunitType , only : lun                
  use PatchType    , only : patch                
  use CNDVType     , only : dgvs_type
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynCNDV_init      ! initialize CNDV weight updates
  public :: dynCNDV_interp    ! interpolate CNDV weight updates to the time step
   !-----------------------------------------------------------------------

contains

   !-----------------------------------------------------------------------
   subroutine dynCNDV_init(bounds, dgvs_inst)
     !
     ! !DESCRIPTION:
     ! Initialize time interpolation of cndv pft weights from annual to time step
     !
     ! Should be called once, in model initialization
     !
     ! !USES:
     use clm_varctl, only : nsrest, nsrStartup
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in)    :: bounds  
     type(dgvs_type)  , intent(inout) :: dgvs_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: ier, p                        ! error status, do-loop index
     character(len=32) :: subname='dynCNDV_init' ! subroutine name
     !-----------------------------------------------------------------------

     if (nsrest == nsrStartup) then
        do p = bounds%begp,bounds%endp
           dgvs_inst%fpcgrid_patch(p) = patch%wtcol(p)
           dgvs_inst%fpcgridold_patch(p) = patch%wtcol(p)
        end do
     end if

  end subroutine dynCNDV_init

  !-----------------------------------------------------------------------
  subroutine dynCNDV_interp( bounds, dgvs_inst)
    !
    ! !DESCRIPTION:
    ! Time interpolate cndv pft weights from annual to time step
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date, get_step_size_real, get_nstep, get_curr_yearfrac
    use landunit_varcon , only : istsoil ! CNDV incompatible with dynLU
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds  
    type(dgvs_type)  , intent(inout) :: dgvs_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,g,l,p            ! indices
    real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: wt1                ! time interpolation weights (weight of time 1)
    real(r8) :: dtime              ! model time step
    real(r8) :: days_per_year      ! days per year
    integer  :: nstep              ! time step number
    integer  :: year               ! year (0, ...) at nstep + 1
    integer  :: mon                ! month (1, ..., 12) at nstep + 1
    integer  :: day                ! day of month (1, ..., 31) at nstep + 1
    integer  :: sec                ! seconds into current date at nstep + 1
    character(len=32) :: subname='dynCNDV_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! Interpolate pft weight to current time step
    ! Map interpolated pctpft to subgrid weights
    ! assumes maxsoil_patches includes bare ground, each landunit has 1 column, 
    ! SCAM not defined and create_croplandunit = .false.

    nstep         = get_nstep()
    dtime         = get_step_size_real()

    wt1 = 1.0_r8 - get_curr_yearfrac(offset = -int(dtime))

    call get_curr_date(year, mon, day, sec, offset=int(dtime))

    do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       l = patch%landunit(p)

       if (lun%itype(l) == istsoil .and. lun%wtgcell(l) > 0._r8) then ! CNDV incompatible with dynLU
          patch%wtcol(p)   = dgvs_inst%fpcgrid_patch(p) + &
                    wt1 * (dgvs_inst%fpcgridold_patch(p) - dgvs_inst%fpcgrid_patch(p))

          if (mon==1 .and. day==1 .and. sec==dtime) then
             dgvs_inst%fpcgridold_patch(p) = dgvs_inst%fpcgrid_patch(p)
          end if
       end if
    end do

  end subroutine dynCNDV_interp

end module dynCNDVMod
