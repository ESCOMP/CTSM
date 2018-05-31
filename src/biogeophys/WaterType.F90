module WaterType

  ! This is the only type that will appear directly in clm_instMod: clm_inst will have a
  ! single instance of water_type, named water_inst, and will call water_inst%Init and
  ! water_inst%Restart.
  type, public :: water_type
     type(waterstate_bulk_type) :: waterstate_bulk_inst
     type(waterdiag_bulk_type) :: waterdiag_bulk_inst
     type(waterbalance_type) :: waterbalance_bulk_inst

     ! Probably move fluxes in here, too - eventually

     ! Mat: You don't need to implement these tracer instances now: These will come in
     ! the next iteration
     type(waterstate_type), allocatable :: waterstate_tracer_inst(:)
     type(waterdiag_type), allocatable :: waterdiag_tracer_inst(:)
     type(waterbalance_type), allocatable :: waterbalance_tracer_inst(:)
   contains
     procedure :: Init
     procedure :: Restart
  end type water_type

contains

  subroutine Init(this, bounds, ...)
    call this%waterstate_bulk_inst%Init(bounds, ...)
    call this%waterdiag_bulk_inst%Init(bounds, ...)
    call this%waterbalance_bulk_inst%Init(bounds, ...)

    ! Mat: You don't need to implement the following now: This will come in the next iteration

    ! First have code to determine number of tracers, then do the following:

    allocate(this%waterstate_tracer_inst(ntracers))
    allocate(this%waterdiag_tracer_inst(ntracers))
    allocate(this%waterbalance_tracer_inst(ntracers))
    do i = 1, size(this%watertracer_inst)
       call this%waterstate_tracer_inst(i)%Init(bounds, ...)
       call this%waterdiag_tracer_inst(i)%Init(bounds, ...)
       call this%waterbalance_tracer_inst(i)%Init(bounds, ...)
    end do
  end subroutine Init

  subroutine Restart(this, bounds, ...)
    call this%waterstate_inst%Restart(bounds, ...)
    call this%waterdiag_inst%Restart(bounds, ...)
    call this%waterbalance_inst%Restart(bounds, ...)

    ! Mat: You don't need to implement the following now: This will come in the next iteration
    do i = 1, size(watertracer_inst)
       call this%waterstate_tracer_inst(i)%Restart(bounds, ...)
       call this%waterdiag_tracer_inst(i)%Restart(bounds, ...)
       call this%waterbalance_tracer_inst(i)%Restart(bounds, ...)
    end do
  end subroutine Restart

end module WaterType
