module ctsm_memcheck

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl, only : iulog
    use spmdMod, only : masterproc
    implicit none
    private

    public :: memcheck
contains

    subroutine memcheck(msg)
        use proc_status_vm, only : prt_vm_status, shr_malloc_trim
        use shr_mem_mod, only : shr_mem_getusage
        use shr_sys_mod, only : shr_sys_flush
        character(len=*), intent(in) :: msg

        real(r8) :: msize, mrss  ! Memory size and resident set size

        call shr_malloc_trim()  ! Make sure the OS trims the memory in response to deallocates

        ! Only output memory on main task as memory usage should be similar between tasks
        if (masterproc) then
           call prt_vm_status('CTSM(Memory check): ' // trim(msg))
           call shr_mem_getusage( msize, mrss, prt=.true.)
           write(iulog,*) '        msize, mrss = ',msize, mrss
           call shr_sys_flush(iulog)
        end if

    end subroutine memcheck

end module ctsm_memcheck