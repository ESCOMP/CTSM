module lilac_utils

  implicit none
  private

  public :: this_clock
  public :: lilac_init_atm2lnd
  public :: lilac_init_lnd2atm
  public :: lilac_atm2lnd
  public :: lilac_lnd2atm

  ! Global index space info for atm data
  ! the HOST ATMOSPHERE is also responsible for filling in the gindex information
  ! this is used to create the distgrid for the mesh in lilac ***
  integer, public, allocatable  :: gindex_atm (:)

  type :: atm2lnd_type
     character(len=128) :: fldname
     real*8, pointer    :: dataptr(:)
     character(len=64)  :: units
     logical            :: provided_by_atm
     logical            :: required_fr_atm 
  end type atm2lnd_type
  type(atm2lnd_type), pointer, public :: atm2lnd(:)

  type :: lnd2atm_type
     character(len=128) :: fldname
     real*8, pointer    :: dataptr(:)
     character(len=64)  :: units
  end type lnd2atm_type
  type(atm2lnd_type), pointer, public :: lnd2atm(:)

  type :: this_clock
     integer, pointer :: yy
     integer, pointer :: mm
     integer, pointer :: dd
     integer, pointer :: hh
     integer, pointer :: mn
     integer, pointer :: ss
  end type this_clock

!========================================================================
contains
!========================================================================

  ! *** NOTE - THE HOST ATMOSPHERE IS RESPONSIBLE for calling
  ! lilac_init that then calls the initialization routines for atm2lnd and lnd2atm

  ! host atm init call will simply be
  ! call lilac_init()

  ! host atm run phase will be
  ! call lilac_atm2lnd(fldname, data1d)

  subroutine lilac_init_atm2lnd(lsize)
    integer, intent(in) :: lsize
    integer :: n

    ! TODO: how is the atm going to specify which fields are not provided = 
    ! should it pass an array of character strings or a colon deliminited set of fields
    ! to specify the fields it will not provide - and then these are checked against those fields

    call atm2lnd_add_fld (atm2lnd, fldname='Sa_z'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Sa_topo'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Sa_u'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Sa_v'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Sa_ptem'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Sa_pbot'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Sa_tbot'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Sa_shum'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_lwdn'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_rainc' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_rainl' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_snowc' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_snowl' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_swndr' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_swvdr' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_swndf' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call atm2lnd_add_fld (atm2lnd, fldname='Faxa_swvdf' , units='unknown', required_fr_atm=.true.  , lsize=lsize)

    ! TODO: optional fields - if these are uncommented then need to make sure that they are also appear in the lnd
    ! import state
    ! CRITICAL the fields in the export state from lilac_atmcap MUST match the fields in the import state to the land 
    ! this is not being checked currently and msut be
    !call atm2lnd_add_fld (atm2lnd, fldname='Sa_methane' , units='unknown', required_fr_atm=.false. , lsize=lsize)
    !call atm2lnd_add_fld (atm2lnd, fldname='Faxa_bcph'  , units='unknown', required_fr_atm=.false. , lsize=lsize)

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

    call lnd2atm_add_fld (lnd2atm, fldname='Sl_lfrin'  , units='unknown', lsize=lsize)
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_t'      , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_tref'   , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_qref'   , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_avsdr'  , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_anidr'  , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_avsdf'  , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_anidf'  , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_snowh'  , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_u10'    , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_fv'     , units='unknown', lsize=lsize) 
    call lnd2atm_add_fld (lnd2atm, fldname='Sl_ram1'   , units='unknown', lsize=lsize) 

    ! TODO: for now are commenting these since they are in the lnd send - however this
    ! is not correct and the lnd send should reintroduce these as soon as possible and
    ! the following should be uncommented
    !call lnd2atm_add_fld (lnd2atm, fldname='Fall_lwup' , units='unknown', lsize=lsize) 
    !call lnd2atm_add_fld (lnd2atm, fldname='Fall_taux' , units='unknown', lsize=lsize) 
    !call lnd2atm_add_fld (lnd2atm, fldname='Fall_tauy' , units='unknown', lsize=lsize) 

    ! now add dataptr memory for all of the fields
    do n = 1,size(lnd2atm)
       allocate(lnd2atm(n)%dataptr(lsize))
    end do
  end subroutine lilac_init_lnd2atm

!========================================================================

  subroutine lilac_atm2lnd(fldname, data)

    ! input/output variables
    character(len=*), intent(in) :: fldname
    real*8, intent(in)           :: data(:)

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
    character(len=*), intent(in) :: fldname
    real*8, intent(out)          :: data(:)

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

  subroutine atm2lnd_add_fld(flds, fldname, units, required_fr_atm, lsize)

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

  end subroutine atm2lnd_add_fld

!========================================================================

  subroutine lnd2atm_add_fld(flds, fldname, units, lsize)

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

  end subroutine lnd2atm_add_fld

end module lilac_utils
