module WaterType

  ! This type encapsulates the management of multiple tracer instances, preventing that
  ! logic from making its way all the way up to clm_inst.
  !
  ! This is the only type that will appear directly in clm_instMod: clm_inst will have a
  ! single instance of water_type, named water_inst, and will call water_inst%Init and
  ! water_inst%Restart.
  type, public :: water_type
     type(waterbulk_type) :: waterbulk_inst

     ! Mat: You don't need to implement this now: This will come in the next iteration
     type(watertracer_type), allocatable :: watertracer_inst(:)
   contains
     procedure :: Init
     procedure :: Restart
  end type water_type

contains

  subroutine Init(this, bounds, ...)
    call this%waterbulk_inst%Init(bounds, ...)

    ! Mat: You don't need to implement the following now: This will come in the next iteration

    ! Determine number of tracers

    allocate(this%watertracer_inst(ntracers))
    do i = 1, size(this%watertracer_inst)
       call this%watertracer_inst(i)%Init(bounds, ...)
    end do
  end subroutine Init

  subroutine Restart(this, bounds, ...)
    call this%waterbulk_inst%Restart(bounds, ...)

    ! Mat: You don't need to implement the following now: This will come in the next iteration
    do i = 1, size(watertracer_inst)
       call this%watertracer_inst(i)%Restart(bounds, ...)
    end do
  end subroutine Restart

end module WaterType
