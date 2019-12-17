module lilac_atmcap

  !-----------------------------------------------------------------------
  ! This is an ESMF lilac cap for the host atmosphere
  !
  ! THE HOST ATMOSPHERE IS RESPONSIBLE for calling lilac_init() and in turn
  ! lilac_init() calls the initialization routines for atm2lnd and lnd2atm
  !
  ! the host atm init call will be
  !      call lilac_init()
  ! the host atm run phase will be
  !     call lilac_atm2lnd(fldname, data1d)
  !     call lilac_run(restart_alarm_is_ringing, stop_alarm_is_ringing)
  !     call lilac_lnd2atm(fldname, data1d)
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod  , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_sys_mod   , only : shr_sys_abort
  use lilac_methods , only : chkerr

  implicit none

  public :: lilac_atmcap_init
  public :: lilac_atmcap_atm2lnd
  public :: lilac_atmcap_lnd2atm
  public :: lilac_atmcap_register

  private :: lilac_atmcap_add_fld

  ! Time invariant input from host atmosphere
  integer, public, allocatable :: gindex_atm(:) ! global index space
  real   , public, allocatable :: atm_lons(:)   ! local longitudes
  real   , public, allocatable :: atm_lats(:)   ! local latitudes
  integer, public              :: atm_global_nx
  integer, public              :: atm_global_ny

  ! Time variant input from host atmosphere
  real(r8) :: nextsw_cday = 1.e36_r8  ! calendar day of the next sw calculation

  type :: atmcap_type
     character(len=CL)  :: fldname
     real(r8), pointer  :: dataptr(:)
     character(len=CS)  :: units
     logical            :: provided_by_atm
     logical            :: required_fr_atm
  end type atmcap_type
  type(atmcap_type), pointer, public :: atm2lnd(:)
  type(atmcap_type), pointer, public :: lnd2atm(:)

  integer                 :: mytask
  integer     , parameter :: debug = 0 ! internal debug level
  character(*), parameter :: u_FILE_u = &
       __FILE__

!========================================================================
contains
!========================================================================

  subroutine lilac_atmcap_init_vars(atm_gindex_in, atm_lons_in, atm_lats_in, atm_global_nx_in, atm_global_ny_in)

    ! input/output variables
    integer , intent(in) :: atm_gindex_in(:)
    real    , intent(in) :: atm_lons_in(:)
    real    , intent(in) :: atm_lats_in(:)
    integer , intent(in) :: atm_global_nx_in
    integer , intent(in) :: atm_global_ny_in

    ! glocal variables
    integer :: n, lsize, fileunit
    ! --------------------------------------------

    lsize = size(atm_gindex_in)
    allocate(gindex_atm(lsize))
    allocate(atm_lons(lsize))
    allocate(atm_lats(lsize))

    ! set module variables
    gindex_atm(:) = atm_gindex_in(:)
    atm_lons(:)   = atm_lons_in(:)
    atm_lats(:)   = atm_lats_in(:)
    atm_global_nx = atm_global_nx_in
    atm_global_ny = atm_global_ny_in

    !-------------------------------------------------------------------------
    ! Set module arrays atm2lnd and lnd2atm
    !-------------------------------------------------------------------------

    ! TODO: how is the atm going to specify which fields are not provided =
    ! should it pass an array of character strings or a colon deliminited set of fields
    ! to specify the fields it will not provide - and then these are checked against those fields

    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_z'          , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_topo'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_u'          , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_v'          , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_ptem'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_pbot'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_tbot'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_shum'       , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_lwdn'     , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_rainc'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_rainl'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_snowc'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_snowl'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_swndr'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_swvdr'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_swndf'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_swvdf'    , units='unknown', required_fr_atm=.true.  , lsize=lsize)

    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_bcphidry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_bcphodry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_bcphiwet' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_ocphidry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_ocphodry' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_ocphiwet' , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstwet1'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstdry1'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstwet2'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstdry2'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstwet3'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstdry3'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstwet4'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
    call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_dstdry4'  , units='unknown', required_fr_atm=.true.  , lsize=lsize)
  ! call lilac_atmcap_add_fld (atm2lnd, fldname='Sa_methane'    , units='unknown', required_fr_atm=.false. , lsize=lsize)
  ! call lilac_atmcap_add_fld (atm2lnd, fldname='Faxa_bcph'     , units='unknown', required_fr_atm=.false. , lsize=lsize)

    ! now add dataptr memory for all of the fields and set default values of provided_by_atm to false
    do n = 1,size(atm2lnd)
       allocate(atm2lnd(n)%dataptr(lsize))
       atm2lnd(n)%provided_by_atm = .false.
    end do

    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_lfrin'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_t'         , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_tref'      , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_qref'      , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_avsdr'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_anidr'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_avsdf'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_anidf'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_snowh'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_u10'       , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_fv'        , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Sl_ram1'      , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_taux'    , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_tauy'    , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_lat'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_sen'     , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_lwup'    , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_evap'    , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_swnet'   , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_flxdst1' , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_flxdst2' , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_flxdst3' , units='unknown', lsize=lsize)
    call lilac_atmcap_add_fld (lnd2atm , fldname='Fall_flxdst4' , units='unknown', lsize=lsize)

    ! now add dataptr memory for all of the fields
    do n = 1,size(lnd2atm)
       allocate(lnd2atm(n)%dataptr(lsize))
    end do

  end subroutine lilac_atmcap_init_vars

  !========================================================================
  subroutine lilac_atmcap_register (comp, rc)

    ! input/output variables
    type(ESMF_GridComp)          :: comp   ! must not be optional
    integer, intent(out)         :: rc

    ! local variables
    type(ESMF_VM)                :: vm
    character(len=*), parameter  :: subname='(lilac_atmcap_register): '
    !-------------------------------------------------------------------------

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
       print *, "in user register routine"
    end if

    ! Initialize return code
    rc = ESMF_SUCCESS

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=lilac_atmcap_init, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=lilac_atmcap_run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=lilac_atmcap_final, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_atmcap_register

!========================================================================

  subroutine lilac_atmcap_init (comp, lnd2atm_state, atm2lnd_state, clock, rc)

    ! input/output variables
    type (ESMF_GridComp) :: comp
    type (ESMF_State)    :: lnd2atm_state
    type (ESMF_State)    :: atm2lnd_state
    type (ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    ! local variables
    integer                :: fileunit
    type(ESMF_Mesh)        :: atm_mesh
    type(ESMF_DistGrid)    :: atm_distgrid
    type(ESMF_Field)       :: field
    type(ESMF_FieldBundle) :: c2a_fb , a2c_fb
    integer                :: n, i, ierr
    integer                :: lsize
    character(len=cl)      :: atm_mesh_filename
    character(len=cl)      :: cvalue
    integer                :: spatialDim
    integer                :: numOwnedElements
    real(r8), pointer      :: ownedElemCoords(:)
    real(r8)               :: mesh_lon, mesh_lat
    real(r8)               :: tolerance = 1.e-5_r8
    character(len=*), parameter :: subname='(lilac_atmcap_init): '
    !-------------------------------------------------------------------------

    namelist /lilac_atmcap_input/ atm_mesh_filename

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"------------------------!", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Read in the atm mesh
    !-------------------------------------------------------------------------

    ! read in mesh file name from namelist
    open(newunit=fileunit, status="old", file="lilac_in")
    read(fileunit, lilac_atmcap_input, iostat=ierr)
    if (ierr > 0) then
       call shr_sys_abort(trim(subname) // 'error reading in lilac_atmcap_input')
    end if
    close(fileunit)

    ! Note that in the call to lilac_atm the host atmospere sent both the gindex_atm and
    ! the atm_mesh_filename that were then set as module variables here

    atm_distgrid = ESMF_DistGridCreate (arbSeqIndexList=gindex_atm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: the addUserArea failed for the 4x5 grid - need to have a
    ! more robust approach - unless the area will simply be ignored for now?
    ! atm_mesh = ESMF_MeshCreate(filename=trim(atm_mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
    !      elementDistGrid=atm_distgrid, addUserArea=.true., rc=rc)
    atm_mesh = ESMF_MeshCreate(filename=trim(atm_mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistGrid=atm_distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) then
       call shr_sys_abort(trim(subname) // 'Error: failure in creating lilac atmcap from meshfile '//&
            trim(atm_mesh_filename))
    end if

    call ESMF_LogWrite(trim(subname)//"Mesh for atmosphere is created for "//trim(atm_mesh_filename), ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // "Mesh for atmosphere is created for "//trim(atm_mesh_filename)
    end if

    !-------------------------------------------------------------------------
    ! Check that lons and lats from the host atmospere match those read
    ! in from the atm mesh file
    !-------------------------------------------------------------------------

    call ESMF_MeshGet(atm_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(gindex_atm)
    if (numOwnedElements /= lsize) then
       print *, 'numOwnedElements in atm_mesh = ',numOwnedElements
       print *, 'local size from gindex_atm from host atm = ',lsize
       call shr_sys_abort('ERROR: numOwnedElements is not equal to lsize')
    end if

    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    call ESMF_MeshGet(atm_mesh, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, lsize
       mesh_lon = ownedElemCoords(2*n-1)
       mesh_lat = ownedElemCoords(2*n)
       if ( abs(mesh_lon - atm_lons(n)) > tolerance) then
          write(6,101),n, atm_lons(n), mesh_lon
101       format('ERROR: lilac_atmcap: n, lon, mesh_lon = ',i6,2(f20.10,2x))
          call shr_sys_abort()
       end if
       if ( abs(mesh_lat - atm_lats(n)) > tolerance) then
          write(6,102),n, atm_lats(n), mesh_lat
102       format('ERROR: lilac_atmcap: n, lat, mesh_lat = ',i6,2(f20.10,2x))
          call shr_sys_abort()
       end if
    end do
    deallocate(ownedElemCoords)

    !-------------------------------------------------------------------------
    ! Create a2c_fb field bundle and add to atm2lnd_state
    !-------------------------------------------------------------------------

    ! create empty field bundle "a2c_fb"
    a2c_fb = ESMF_FieldBundleCreate(name="a2c_fb", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create fields and add to field bundle
    do n = 1, size(atm2lnd)
       field = ESMF_FieldCreate(atm_mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
            name=trim(atm2lnd(n)%fldname), farrayPtr=atm2lnd(n)%dataptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleAdd(a2c_fb, (/field/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (debug > 0) then
          call ESMF_FieldPrint(field,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    ! add field bundle to atm2lnd_state
    call ESMF_StateAdd(atm2lnd_state, (/a2c_fb/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"lilac a2c_fb fieldbundle created and added to atm2lnd_state", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "lilac a2c_fb fieldbundle created and added to atm2lnd_state"
    end if

    ! add nextsw_cday attributes
    write(cvalue,*) nextsw_cday
    call ESMF_AttributeSet(atm2lnd_state, name="nextsw_cday", value=cvalue, rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------------------------------------
    ! Create c2a_fb field bundle and add to lnd2atm_state
    !-------------------------------------------------------------------------

    ! create empty field bundle "c2a_fb"
    c2a_fb = ESMF_FieldBundleCreate (name="c2a_fb", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create fields and add to field bundle
    do n = 1, size(lnd2atm)
       field = ESMF_FieldCreate(atm_mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
            name=trim(lnd2atm(n)%fldname), farrayPtr=lnd2atm(n)%dataptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleAdd(c2a_fb, (/field/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (debug > 0) then
          call ESMF_FieldPrint(field,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    ! add field bundle to lnd2atm_state
    call ESMF_StateAdd(lnd2atm_state, (/c2a_fb/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"lilac c2a_fb fieldbundle is done and added to lnd2atm_state", ESMF_LOGMSG_INFO)

  end subroutine lilac_atmcap_init

!========================================================================
  subroutine lilac_atmcap_run(comp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Initialize return code
    ! This routine does nothing - its all in the atm->lnd coupler
    rc = ESMF_SUCCESS

  end subroutine lilac_atmcap_run

!========================================================================
  subroutine lilac_atmcap_final(comp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle) ::  import_fieldbundle, export_fieldbundle
    character(len=*), parameter :: subname='( lilac_atmcap_final): '
    !-------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_StateGet(importState, "c2a_fb", import_fieldbundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(exportState, "a2c_fb", export_fieldbundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleDestroy(import_fieldbundle, rc=rc)
    call ESMF_FieldBundleDestroy(export_fieldbundle, rc=rc)

    call ESMF_LogWrite(subname//"Finished lilac_atmcap_final", ESMF_LOGMSG_INFO)

  end subroutine lilac_atmcap_final

!========================================================================
  subroutine lilac_atmcap_atm2lnd(fldname, data)

    ! input/output variables
    character(len=*), intent(in) :: fldname
    real(r8), intent(in)         :: data(:)

    ! local variables
    integer :: n
    logical :: found
    character(len=*), parameter  :: subname='(lilac_atmcap_atm2lnd)'
    ! --------------------------------------------

    found = .false.
    do n = 1,size(atm2lnd)
       if (trim(fldname) == atm2lnd(n)%fldname) then
          found = .true.
          if (size(data) /= size(atm2lnd(n)%dataptr)) then
             call shr_sys_abort(trim(subname) // 'size(data) not equal to size(atm2lnd(n)%dataptr')
          else
             atm2lnd(n)%dataptr(:) = data(:)
          end if
          atm2lnd(n)%provided_by_atm = .true.
          exit
       end if
    end do
    if (.not. found) then
       call shr_sys_abort(trim(subname) // 'atm2lnd field name ' // trim(fldname) //' not found')
    end if

  contains

    subroutine lilac_atm2lnd_check()
      integer :: n    ! if there are fields that the atmosphere does not provide but
      ! that are required - then abort
      do n = 1,size(atm2lnd)
         if (atm2lnd(n)%required_fr_atm .and. (.not. atm2lnd(n)%provided_by_atm)) then
            ! call abort or provide default values?
         else if (.not. atm2lnd(n)%provided_by_atm) then
            ! create default values
         end if
      end do
    end subroutine lilac_atm2lnd_check

  end subroutine lilac_atmcap_atm2lnd

!========================================================================
  subroutine lilac_atmcap_lnd2atm(fldname, data)

    ! input/output variables
    character(len=*) , intent(in)  :: fldname
    real(r8)         , intent(out) :: data(:)

    ! local variables
    integer :: n
    character(len=*), parameter  :: subname='(lilac_atmcap_lnd2atm)'
    ! --------------------------------------------

    do n = 1,size(lnd2atm)
       if (trim(fldname) == lnd2atm(n)%fldname) then
          if (size(data) /= size(lnd2atm(n)%dataptr)) then
             call shr_sys_abort(trim(subname) // 'size(data) not equal to size(lnd2atm(n)%dataptr')
          else
             data(:) = lnd2atm(n)%dataptr(:)
          end if
       end if
    end do
  end subroutine lilac_atmcap_lnd2atm

!========================================================================
  subroutine lilac_atmcap_add_fld(flds, fldname, units, lsize, required_fr_atm)

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

    type(atmcap_type), pointer   :: flds(:)
    character(len=*) , intent(in) :: fldname
    character(len=*) , intent(in) :: units
    integer          , intent(in) :: lsize
    logical, optional, intent(in) :: required_fr_atm

    ! local variables
    integer :: n,oldsize,newsize
    type(atmcap_type), pointer :: newflds(:)
    character(len=*), parameter :: subname='(lilac_atmcap_atm2lnd_fld)'
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
       if (present(required_fr_atm)) then
          flds(newsize)%required_fr_atm = required_fr_atm
       end if

    else
       allocate(flds(newsize))
       flds(newsize)%fldname   = trim(fldname)
       flds(newsize)%units     = trim(units)
       if (present(required_fr_atm)) then
          flds(newsize)%required_fr_atm = required_fr_atm
       end if
    end if

  end subroutine lilac_atmcap_add_fld

end module lilac_atmcap
