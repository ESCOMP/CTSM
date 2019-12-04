module lilac_utils

  ! ***********************************************************************
  ! THE HOST ATMOSPHERE IS RESPONSIBLE for calling lilac_init() and in turn
  ! lilac_init() calls the initialization routines for atm2lnd and lnd2atm
  !
  ! the host atm init call will be 
  !      call lilac_init()
  ! the host atm run phase will be 
  !     call lilac_atm2lnd(fldname, data1d)
  !     call lilac_run(restart_alarm_is_ringing, stop_alarm_is_ringing)
  !     call lilac_lnd2atm(fldname, data1d)
  ! ***********************************************************************

  use ESMF
  use shr_kind_mod  , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_sys_mod   , only : shr_sys_abort
  use lilac_methods , only : chkerr

  implicit none
  private

  public :: lilac_init_atm2lnd
  public :: lilac_init_lnd2atm
  public :: lilac_atm2lnd
  public :: lilac_lnd2atm
  public :: lilac_field_bundle_to_land
  public :: lilac_field_bundle_fr_land

  private :: lilac_atm2lnd_add_fld
  private :: lilac_lnd2atm_add_fld

  ! Global index space info for atm data
  integer, public, allocatable  :: gindex_atm (:)

  ! Mesh file to be read in by lilac_atm
  character(len=CL), public :: atm_mesh_filename

  ! Mesh file to be read in by ctsm
  character(len=CL), public :: lnd_mesh_filename

  type :: atm2lnd_type
     character(len=CL) :: fldname
     real(r8), pointer  :: dataptr(:)
     character(len=CS)  :: units
     logical            :: provided_by_atm
     logical            :: required_fr_atm
  end type atm2lnd_type
  type(atm2lnd_type), pointer, public :: atm2lnd(:)

  type :: lnd2atm_type
     character(len=128) :: fldname
     real(r8), pointer  :: dataptr(:)
     character(len=64)  :: units
  end type lnd2atm_type
  type(atm2lnd_type), pointer, public :: lnd2atm(:)

  character(*), parameter :: u_FILE_u = &
       __FILE__

!========================================================================
contains
!========================================================================

  subroutine lilac_init_atm2lnd(lsize)
    integer, intent(in) :: lsize
    integer :: n

    ! TODO: how is the atm going to specify which fields are not provided =
    ! should it pass an array of character strings or a colon deliminited set of fields
    ! to specify the fields it will not provide - and then these are checked against those fields

    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_z'          , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_topo'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_u'          , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_v'          , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_ptem'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_pbot'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_tbot'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_shum'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_lwdn'     , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_rainc'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_rainl'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_snowc'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_snowl'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_swndr'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_swvdr'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_swndf'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_swvdf'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)

    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_bcphidry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_bcphodry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_bcphiwet' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_ocphidry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_ocphodry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_ocphiwet' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstwet1'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstdry1'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstwet2'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstdry2'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstwet3'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstdry3'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstwet4'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_dstdry4'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)

    ! TODO: optional fields - if these are uncommented then need to make sure that they are also appear in the lnd
    ! import state
    ! CRITICAL the fields in the export state from lilac_atmcap MUST match the fields in the import state to the land
    ! this is not being checked currently and msut be
    !call lilac_atm2lnd_add_fld (atm2lnd, fldname='Sa_methane' , units='unknown', required_fr_atm=.false. , lsize=lsize)
    !call lilac_atm2lnd_add_fld (atm2lnd, fldname='Faxa_bcph'  , units='unknown', required_fr_atm=.false. , lsize=lsize)

    ! now add dataptr memory for all of the fields and set default values of provided_by_atm to false
    do n = 1,size(atm2lnd)
       allocate(atm2lnd(n)%dataptr(lsize))
       atm2lnd(n)%provided_by_atm = .false.
    end do
  end subroutine lilac_init_atm2lnd

  !========================================================================

  subroutine lilac_init_lnd2atm(lsize)
    integer, intent(in) :: lsize
    integer :: n

    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_lfrin'  , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_t'      , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_tref'   , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_qref'   , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_avsdr'  , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_anidr'  , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_avsdf'  , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_anidf'  , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_snowh'  , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_u10'    , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_fv'     , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Sl_ram1'   , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Fall_lwup' , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Fall_taux' , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Fall_tauy' , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Fall_evap' , units='unknown', lsize=lsize)
    call lilac_lnd2atm_add_fld (lnd2atm, fldname='Fall_swnet', units='unknown', lsize=lsize)

    ! now add dataptr memory for all of the fields
    do n = 1,size(lnd2atm)
       allocate(lnd2atm(n)%dataptr(lsize))
    end do
  end subroutine lilac_init_lnd2atm

!========================================================================

  subroutine lilac_atm2lnd(fldname, data)

    ! input/output variables
    character(len=*), intent(in) :: fldname
    real(r8), intent(in)         :: data(:)

    ! local variables
    integer :: n
    logical :: found
    ! --------------------------------------------

    found = .false.
    do n = 1,size(atm2lnd)
       if (trim(fldname) == atm2lnd(n)%fldname) then
          found = .true.
          if (size(data) /= size(atm2lnd(n)%dataptr)) then
             ! call abort - TODO: what is the abort call in lilac
          else
             atm2lnd(n)%dataptr(:) = data(:)
          end if
          atm2lnd(n)%provided_by_atm = .true.
          exit
       end if
    end do
    if (.not. found) then
       ! abort
    end if

  end subroutine lilac_atm2lnd

!========================================================================

  subroutine lilac_atm2lnd_check()

    ! local variables
    integer :: n
    ! --------------------------------------------

    ! if there are fields that the atmosphere does not provide but that are required - then abort
    do n = 1,size(atm2lnd)
       if (atm2lnd(n)%required_fr_atm .and. (.not. atm2lnd(n)%provided_by_atm)) then
          ! call abort or provide default values?
       else if (.not. atm2lnd(n)%provided_by_atm) then
          ! create default values
       end if
    end do
  end subroutine lilac_atm2lnd_check

!========================================================================

  subroutine lilac_lnd2atm(fldname, data)
    ! input/output variables
    character(len=*) , intent(in)  :: fldname
    real(r8)         , intent(out) :: data(:)

    ! local variables
    integer :: n
    ! --------------------------------------------

    do n = 1,size(lnd2atm)
       if (trim(fldname) == lnd2atm(n)%fldname) then
          if (size(data) /= size(lnd2atm(n)%dataptr)) then
             ! call abort - TODO: what is the abort call in lilac
          else
             data(:) = lnd2atm(n)%dataptr(:)
          end if
       end if
    end do
  end subroutine lilac_lnd2atm

!========================================================================

  subroutine lilac_atm2lnd_add_fld(flds, fldname, units, required_fr_atm, lsize)

    ! ----------------------------------------------
    ! Add an entry to to the flds array
    ! Use pointers to create an extensible allocatable array.
    ! to allow the size of flds to grow, the process for
    ! adding a new field is:
    ! 1) allocate newflds to be N (one element larger than flds)
    ! 2) copy flds into first N-1 elements of newflds
    ! 3) newest flds entry is Nth element of newflds
    ! 4) deallocate / nullify flds
    ! 5) point flds => newflds
    ! ----------------------------------------------

    type(atm2lnd_type), pointer   :: flds(:)
    character(len=*) , intent(in) :: fldname
    character(len=*) , intent(in) :: units
    logical          , intent(in) :: required_fr_atm
    integer          , intent(in) :: lsize

    ! local variables
    integer :: n,oldsize,newsize
    type(atm2lnd_type), pointer :: newflds(:)
    character(len=*), parameter :: subname='(lilac_utils_add_atm2lnd_fld)'
    ! ----------------------------------------------

    if (associated(flds)) then
       oldsize = size(flds)
    else
       oldsize = 0
    end if
    newsize = oldsize + 1

    if (oldsize > 0) then
       ! 1) allocate newfld to be size (one element larger than input flds)
       allocate(newflds(newsize))

       ! 2) copy flds into first N-1 elements of newflds
       do n = 1,oldsize
          newflds(n)%fldname    =  flds(n)%fldname
          newflds(n)%units      =  flds(n)%units
          newflds(n)%required_fr_atm = flds(n)%required_fr_atm
       end do

       ! 3) deallocate / nullify flds
       if (oldsize >  0) then
          deallocate(flds)
          nullify(flds)
       end if

       ! 4) point flds => new_flds
       flds => newflds

       ! 5) update flds information for new entry
       flds(newsize)%fldname   = trim(fldname)
       flds(newsize)%units     = trim(units)
       flds(newsize)%required_fr_atm = required_fr_atm

    else
       allocate(flds(newsize))
       flds(newsize)%fldname   = trim(fldname)
       flds(newsize)%units     = trim(units)
       flds(newsize)%required_fr_atm = required_fr_atm
    end if

  end subroutine lilac_atm2lnd_add_fld

!========================================================================

  subroutine lilac_lnd2atm_add_fld(flds, fldname, units, lsize)

    ! ----------------------------------------------
    ! Add an entry to to the flds array
    ! Use pointers to create an extensible allocatable array.
    ! to allow the size of flds to grow, the process for
    ! adding a new field is:
    ! 1) allocate newflds to be N (one element larger than flds)
    ! 2) copy flds into first N-1 elements of newflds
    ! 3) newest flds entry is Nth element of newflds
    ! 4) deallocate / nullify flds
    ! 5) point flds => newflds
    ! ----------------------------------------------

    type(atm2lnd_type), pointer   :: flds(:)
    character(len=*) , intent(in) :: fldname
    character(len=*) , intent(in) :: units
    integer          , intent(in) :: lsize

    ! local variables
    integer :: n,oldsize,newsize
    type(atm2lnd_type), pointer :: newflds(:)
    character(len=*), parameter :: subname='(lilac_init_lnd2atm)'
    ! ----------------------------------------------

    if (associated(flds)) then
       oldsize = size(flds)
    else
       oldsize = 0
    end if
    newsize = oldsize + 1

    ! 1) allocate newfld to be size (one element larger than input flds)
    allocate(newflds(newsize))

    ! 2) copy flds into first N-1 elements of newflds
    do n = 1,oldsize
       newflds(n)%fldname    = flds(n)%fldname
       newflds(n)%units      = flds(n)%units
    end do

    ! 3) deallocate / nullify flds
    if (oldsize >  0) then
       deallocate(flds)
       nullify(flds)
    end if

    ! 4) point flds => new_flds
    flds => newflds

    ! 5) now update flds information for new entry
    flds(newsize)%fldname   = trim(fldname)
    flds(newsize)%units     = trim(units)

  end subroutine lilac_lnd2atm_add_fld

!========================================================================

  subroutine lilac_field_bundle_to_land(mesh, fieldbundle, rc)

    type(ESMF_Mesh)        :: mesh
    type(ESMF_FieldBundle) :: fieldbundle
    integer,  intent(out)  :: rc

    integer :: n

    rc = ESMF_SUCCESS

    ! Add empty fields to field bundle
    do n = 1, size(atm2lnd)
       call fldbundle_add(trim(atm2lnd(n)%fldname), mesh, fieldbundle, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine lilac_field_bundle_to_land

!===============================================================================

  subroutine lilac_field_bundle_fr_land(mesh, fieldbundle, rc)

    type(ESMF_Mesh)        :: mesh 
    type(ESMF_FieldBundle) :: fieldbundle
    integer, intent(out)   :: rc

    integer :: n

    rc = ESMF_SUCCESS

    ! Add empty fields to field bundle
    do n = 1, size(lnd2atm)
       call fldbundle_add( trim(lnd2atm(n)%fldname), mesh, fieldbundle, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine lilac_field_bundle_fr_land

!===============================================================================

  subroutine fldbundle_add(fldname, mesh, fieldbundle, rc)

    !---------------------------
    ! Create an empty input field with name 'stdname' to add to fieldbundle
    !---------------------------

    ! input/output variables
    character(len=*)       , intent(in)    :: fldname
    type(ESMF_Mesh)        , intent(in)    :: mesh
    type(ESMF_FieldBundle) , intent(inout) :: fieldbundle
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field) :: field
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldname), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleAdd(fieldbundle, (/field/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine fldbundle_add

end module lilac_utils
