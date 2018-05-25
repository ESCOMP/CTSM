module WaterType

  ! This type encapsulates the management of multiple tracer instances, preventing that
  ! logic from making its way all the way up to clm_inst.
  !
  ! This is the only type that will appear directly in clm_instMod: clm_inst will have a
  ! single instance of water_type, named water_inst, and will call water_inst%Init and
  ! water_inst%Restart.
  type, public :: water_type
     type(water_bulk_and_tracer_type) :: water_bulk_inst
     type(water_bulk_and_tracer_type), allocatable :: water_tracer_inst(:)
     type(water_bulk_only_type) :: water_bulk_only_inst
   contains
     procedure :: Init
     procedure :: Restart
  end type water_type

contains

  subroutine Init(this, bounds, ...)
    call this%water_bulk_inst%Init(bounds, ...)
    call this%water_bulk_only_inst%Init(bounds, ...)

    ! Determine number of tracers

    allocate(this%water_tracer_inst(ntracers))
    do i = 1, size(this%water_tracer_inst)
       call this%water_tracer_inst(i)%Init(bounds, ...)
    end do
  end subroutine Init

  subroutine Restart(this, bounds, ...)
    call this%water_bulk_inst%Restart(bounds, ...)
    call this%water_bulk_only_inst%Restart(bounds, ...)

    do i = 1, size(water_tracer_inst)
       call this%water_tracer_inst(i)%Restart(bounds, ...)
    end do
  end subroutine Restart

end module WaterType
