module dynEDMod
  !
  ! !USES:
  use clmtype
  use decompMod,  only : bounds_type
  use landunit_varcon, only : istsoil
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

  public :: dyn_ED     ! transfers weights calculated internally by ED into wtcol. 
 
contains

  subroutine dyn_ED( bounds )

     use EDclmType , only : EDpft

     type(bounds_type), intent(in) :: bounds  ! bounds
    
     ! !LOCAL VARIABLES:
     integer  ::  p,c           ! indices
  
     associate(&
        pcolumn  => pft%column   &
     )

     do p = bounds%begp,bounds%endp
        c = pcolumn(p)
        if (col%itype(c) == istsoil) then 
          if ((EDpft%ED_patch(p) == 1 ) .or. (EDpft%ED_bareground(p) == 1)) then
            pft%wtcol(p) = EDpft%wtED(p)
          else
            pft%wtcol(p)  = 0.0_r8 
          end if
        end if

     end do
    
  end associate

  end subroutine dyn_ED

end module dynEDMod
