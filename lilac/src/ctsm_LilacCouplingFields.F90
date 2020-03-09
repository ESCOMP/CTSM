module ctsm_LilacCouplingFields

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines the coupling fields between atmosphere and land
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use ctsm_LilacCouplingFieldIndices
  use ctsm_LilacLnd2AtmFieldListType, only : lilac_lnd2atm_field_list_type
  use ctsm_LilacAtm2LndFieldListType, only : lilac_atm2lnd_field_list_type

  implicit none
  private

  !
  ! !PUBLIC ROUTINES:

  ! ------------------------------------------------------------------------
  ! Routines that should be called by the host atmosphere to set / get coupling fields
  ! ------------------------------------------------------------------------

  public :: lilac_atm2lnd  ! Set a single atm -> lnd field
  public :: lilac_lnd2atm  ! Get a single lnd -> atm field

  ! ------------------------------------------------------------------------
  ! Routines that should be used internally by LILAC, *not* called directly from the host
  ! atmosphere
  ! ------------------------------------------------------------------------

  public :: create_a2l_field_list
  public :: create_l2a_field_list
  public :: complete_a2l_field_list
  public :: complete_l2a_field_list

  !
  ! !PUBLIC DATA:

  ! ------------------------------------------------------------------------
  ! These variables should only be used internally by LILAC. The host atmosphere model
  ! should interact with them via the lilac_atm2lnd and lilac_lnd2atm routines.
  ! ------------------------------------------------------------------------

  type(lilac_atm2lnd_field_list_type), public :: a2l_fields
  type(lilac_lnd2atm_field_list_type), public :: l2a_fields

contains

  !-----------------------------------------------------------------------
  subroutine lilac_atm2lnd(field_index, data)
    !
    ! !DESCRIPTION:
    ! Set a single atm -> lnd field
    !
    ! field_index should be one of the lilac_a2l_* indices defined in ctsm_LilacCouplingFieldIndices
    !
    ! !ARGUMENTS:
    integer, intent(in) :: field_index
    real(r8), intent(in) :: data(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'lilac_atm2lnd'
    !-----------------------------------------------------------------------

    call a2l_fields%set_field(field_index, data)

  end subroutine lilac_atm2lnd

  !-----------------------------------------------------------------------
  subroutine lilac_lnd2atm(field_index, data)
    !
    ! !DESCRIPTION:
    ! Get a single lnd -> atm field
    !
    ! field_index should be one of the lilac_l2a_* indices defined in ctsm_LilacCouplingFieldIndices
    !
    ! !ARGUMENTS:
    integer, intent(in) :: field_index
    real(r8), intent(out) :: data(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'lilac_lnd2atm'
    !-----------------------------------------------------------------------

    call l2a_fields%get_field(field_index, data)

  end subroutine lilac_lnd2atm

  !-----------------------------------------------------------------------
  subroutine create_a2l_field_list()
    !
    ! !DESCRIPTION:
    ! Create the list of fields passed from atm -> lnd.
    !
    ! All of the lilac_a2l_* indices are valid after this is called. However, note that
    ! a2l_fields still isn't fully usable until complete_a2l_field_list is called.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'create_a2l_field_list'
    !-----------------------------------------------------------------------

    call a2l_fields%init()

    call a2l_fields%add_var(fieldname='Sa_landfrac'   , units='fraction', available_from_data=.false., &
         can_be_time_const=.true., field_index=lilac_a2l_Sa_landfrac)
    call a2l_fields%add_var(fieldname='Sa_z'          , units='unknown', available_from_data=.false., &
         can_be_time_const=.true., field_index=lilac_a2l_Sa_z)
    call a2l_fields%add_var(fieldname='Sa_topo'       , units='unknown', available_from_data=.false., &
         can_be_time_const=.true., field_index=lilac_a2l_Sa_topo)
    call a2l_fields%add_var(fieldname='Sa_u'          , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Sa_u)
    call a2l_fields%add_var(fieldname='Sa_v'          , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Sa_v)
    call a2l_fields%add_var(fieldname='Sa_ptem'       , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Sa_ptem)
    call a2l_fields%add_var(fieldname='Sa_pbot'       , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Sa_pbot)
    call a2l_fields%add_var(fieldname='Sa_tbot'       , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Sa_tbot)
    call a2l_fields%add_var(fieldname='Sa_shum'       , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Sa_shum)
    call a2l_fields%add_var(fieldname='Faxa_lwdn'     , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_lwdn)
    call a2l_fields%add_var(fieldname='Faxa_rainc'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_rainc)
    call a2l_fields%add_var(fieldname='Faxa_rainl'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_rainl)
    call a2l_fields%add_var(fieldname='Faxa_snowc'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_snowc)
    call a2l_fields%add_var(fieldname='Faxa_snowl'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_snowl)
    call a2l_fields%add_var(fieldname='Faxa_swndr'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_swndr)
    call a2l_fields%add_var(fieldname='Faxa_swvdr'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_swvdr)
    call a2l_fields%add_var(fieldname='Faxa_swndf'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_swndf)
    call a2l_fields%add_var(fieldname='Faxa_swvdf'    , units='unknown', available_from_data=.false., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_swvdf)

    call a2l_fields%add_var(fieldname='Faxa_bcphidry' , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_bcphidry)
    call a2l_fields%add_var(fieldname='Faxa_bcphodry' , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_bcphodry)
    call a2l_fields%add_var(fieldname='Faxa_bcphiwet' , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_bcphiwet)
    call a2l_fields%add_var(fieldname='Faxa_ocphidry' , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_ocphidry)
    call a2l_fields%add_var(fieldname='Faxa_ocphodry' , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_ocphodry)
    call a2l_fields%add_var(fieldname='Faxa_ocphiwet' , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_ocphiwet)
    call a2l_fields%add_var(fieldname='Faxa_dstwet1'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstwet1)
    call a2l_fields%add_var(fieldname='Faxa_dstdry1'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstdry1)
    call a2l_fields%add_var(fieldname='Faxa_dstwet2'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstwet2)
    call a2l_fields%add_var(fieldname='Faxa_dstdry2'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstdry2)
    call a2l_fields%add_var(fieldname='Faxa_dstwet3'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstwet3)
    call a2l_fields%add_var(fieldname='Faxa_dstdry3'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstdry3)
    call a2l_fields%add_var(fieldname='Faxa_dstwet4'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstwet4)
    call a2l_fields%add_var(fieldname='Faxa_dstdry4'  , units='unknown', available_from_data=.true., &
         can_be_time_const=.false., field_index=lilac_a2l_Faxa_dstdry4)

  end subroutine create_a2l_field_list

  !-----------------------------------------------------------------------
  subroutine create_l2a_field_list()
    !
    ! !DESCRIPTION:
    ! Create the list of fields passed from lnd -> atm.
    !
    ! All of the lilac_l2a_* indices are valid after this is called. However, note that
    ! l2a_fields still isn't fully usable until complete_l2a_field_list is called.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'create_l2a_field_list'
    !-----------------------------------------------------------------------

    call l2a_fields%init()

    call l2a_fields%add_var(fieldname='Sl_t'         , units='unknown', &
         field_index=lilac_l2a_Sl_t)
    call l2a_fields%add_var(fieldname='Sl_tref'      , units='unknown', &
         field_index=lilac_l2a_Sl_tref)
    call l2a_fields%add_var(fieldname='Sl_qref'      , units='unknown', &
         field_index=lilac_l2a_Sl_qref)
    call l2a_fields%add_var(fieldname='Sl_avsdr'     , units='unknown', &
         field_index=lilac_l2a_Sl_avsdr)
    call l2a_fields%add_var(fieldname='Sl_anidr'     , units='unknown', &
         field_index=lilac_l2a_Sl_anidr)
    call l2a_fields%add_var(fieldname='Sl_avsdf'     , units='unknown', &
         field_index=lilac_l2a_Sl_avsdf)
    call l2a_fields%add_var(fieldname='Sl_anidf'     , units='unknown', &
         field_index=lilac_l2a_Sl_anidf)
    call l2a_fields%add_var(fieldname='Sl_snowh'     , units='unknown', &
         field_index=lilac_l2a_Sl_snowh)
    call l2a_fields%add_var(fieldname='Sl_u10'       , units='unknown', &
         field_index=lilac_l2a_Sl_u10)
    call l2a_fields%add_var(fieldname='Sl_fv'        , units='unknown', &
         field_index=lilac_l2a_Sl_fv)
    call l2a_fields%add_var(fieldname='Sl_ram1'      , units='unknown', &
         field_index=lilac_l2a_Sl_ram1)
    call l2a_fields%add_var(fieldname='Sl_z0m'       , units='m'      , &
         field_index=lilac_l2a_Sl_z0m)
    call l2a_fields%add_var(fieldname='Fall_taux'    , units='unknown', &
         field_index=lilac_l2a_Fall_taux)
    call l2a_fields%add_var(fieldname='Fall_tauy'    , units='unknown', &
         field_index=lilac_l2a_Fall_tauy)
    call l2a_fields%add_var(fieldname='Fall_lat'     , units='unknown', &
         field_index=lilac_l2a_Fall_lat)
    call l2a_fields%add_var(fieldname='Fall_sen'     , units='unknown', &
         field_index=lilac_l2a_Fall_sen)
    call l2a_fields%add_var(fieldname='Fall_lwup'    , units='unknown', &
         field_index=lilac_l2a_Fall_lwup)
    call l2a_fields%add_var(fieldname='Fall_evap'    , units='unknown', &
         field_index=lilac_l2a_Fall_evap)
    call l2a_fields%add_var(fieldname='Fall_swnet'   , units='unknown', &
         field_index=lilac_l2a_Fall_swnet)
    call l2a_fields%add_var(fieldname='Fall_flxdst1' , units='unknown', &
         field_index=lilac_l2a_Fall_flxdst1)
    call l2a_fields%add_var(fieldname='Fall_flxdst2' , units='unknown', &
         field_index=lilac_l2a_Fall_flxdst2)
    call l2a_fields%add_var(fieldname='Fall_flxdst3' , units='unknown', &
         field_index=lilac_l2a_Fall_flxdst3)
    call l2a_fields%add_var(fieldname='Fall_flxdst4' , units='unknown', &
         field_index=lilac_l2a_Fall_flxdst4)

  end subroutine create_l2a_field_list

  !-----------------------------------------------------------------------
  subroutine complete_a2l_field_list(lsize_atm, fields_needed_from_data)
    !
    ! !DESCRIPTION:
    ! Complete the setup of a2l_fields.
    !
    ! This is separated from create_a2l_field_list because lsize may not be available at
    ! the point when that routine is called. Also, note that this sets
    ! fields_needed_from_data, which won't be available until later in initialization.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: lsize_atm  ! number of atm points on this proc

    ! List of field indices that need to be read from data, because the host atmosphere
    ! isn't going to provide them. These should be indices given in
    ! ctsm_LilacCouplingFields (lilac_a2l_Faxa_bcphidry). This can be an empty list if no
    ! fields need to be read from data.
    integer          , intent(in)    :: fields_needed_from_data(:)

    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'complete_a2l_field_list'
    !-----------------------------------------------------------------------

    call a2l_fields%complete_setup(lsize_atm)
    call a2l_fields%set_needed_from_data(fields_needed_from_data)

  end subroutine complete_a2l_field_list

  !-----------------------------------------------------------------------
  subroutine complete_l2a_field_list(lsize_atm)
    !
    ! !DESCRIPTION:
    ! Complete the setup of l2a_fields.
    !
    ! This is separated from create_l2a_field_list because lsize may not be available at
    ! the point when that routine is called.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: lsize_atm  ! number of atm points on this proc
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'complete_l2a_field_list'
    !-----------------------------------------------------------------------

    call l2a_fields%complete_setup(lsize_atm)

  end subroutine complete_l2a_field_list

end module ctsm_LilacCouplingFields
