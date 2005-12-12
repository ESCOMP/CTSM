#include <misc.h>
#include <preproc.h>

module lnd2atmMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: lnd2atmMod
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd2atm
!
! !INTERFACE: subroutine lnd2atm(init)
  subroutine lnd2atm(init)
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
    use subgridAveMod
    use decompMod   , only : get_proc_bounds
    use clm_varcon  , only : sb
    use clm_varpar  , only : numrad
!
! !ARGUMENTS:
    implicit none
    logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
!
! !REVISION HISTORY:
! Mariana Vertenstein: created 03/10-25
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: begp, endp      ! per-proc beginning and ending pft indices
    integer :: begc, endc      ! per-proc beginning and ending column indices
    integer :: begl, endl      ! per-proc beginning and ending landunit indices
    integer :: begg, endg      ! per-proc gridcell ending gridcell indices
!
! !USES:
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
! 03-08-25 : Updated to vector data structure (Mariana Vertenstein)
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g                          ! indices
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine processor bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Compute gridcell averages. 

   if (present(init)) then

      call c2g(begc, endc, begl, endl, begg, endg, &
           cptr%cws%h2osno, gptr%l2as%h2osno,&
           c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
      do g = begg,endg
         gptr%l2as%h2osno(g) = gptr%l2as%h2osno(g)/1000._r8
      end do
      
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
           pptr%pps%albd, gptr%l2as%albd,&
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
      
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
           pptr%pps%albi, gptr%l2as%albi,&
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
      
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pef%eflx_lwrad_out, gptr%l2af%eflx_lwrad_out,&
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
      do g = begg,endg
         gptr%l2as%t_rad(g) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb))
      end do

   else

      call c2g(begc, endc, begl, endl, begg, endg, cptr%cws%h2osno, gptr%l2as%h2osno,&
           c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
      do g = begg,endg
         gptr%l2as%h2osno(g) = gptr%l2as%h2osno(g)/1000._r8
      end do

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
           pptr%pps%albd, gptr%l2as%albd, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
           pptr%pps%albi, gptr%l2as%albi, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pes%t_ref2m, gptr%l2as%t_ref2m, & 
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pes%q_ref2m, gptr%l2as%q_ref2m, & 
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pmf%taux, gptr%l2af%taux, & 
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pmf%tauy, gptr%l2af%tauy, & 
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pef%eflx_lh_tot, gptr%l2af%eflx_lh_tot, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pef%eflx_sh_tot, gptr%l2af%eflx_sh_tot, & 
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pwf%qflx_evap_tot, gptr%l2af%qflx_evap_tot, & 
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pef%fsa, gptr%l2af%fsa, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
                  
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pef%eflx_lwrad_out, gptr%l2af%eflx_lwrad_out, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
                  
#if (defined CN)
      call c2g(begc, endc, begl, endl, begg, endg, &
           cptr%ccf%nee, gptr%l2af%nee, &
           c2l_scale_type= 'unity', l2g_scale_type='unity')
#elif (defined CASA)
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pps%co2flux, gptr%l2af%nee, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#else
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pcf%fco2, gptr%l2af%nee, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
       ! Note that fco2 in is umolC/m2/sec so units need to be changed to gC/m2/sec
       do g = begg,endg
          gptr%l2af%nee(g) = gptr%l2af%nee(g)*12.011e-6_r8
       end do
#endif
       ! Convert from gC/m2/s to kgC/m2/s
!dir$ concurrent
!cdir nodep
       do g = begg,endg
          gptr%l2af%nee(g) = gptr%l2af%nee(g)*1.0e-3_r8
       end do

!dir$ concurrent
!cdir nodep
      do g = begg,endg
         gptr%l2as%t_rad(g) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb))
      end do

   end if

 end subroutine lnd2atm

end module lnd2atmMod
