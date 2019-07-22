module mksoilUtilsMod

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: mksoilUtils
  !
  ! !DESCRIPTION:
  ! Lower-level utilities used in making soil data.
  !
  ! These are separated out from mksoilMod mainly as an aid to testing.
  !
  ! !REVISION HISTORY:
  ! Author: Bill Sacks
  !
  !-----------------------------------------------------------------------
  !!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkgridmapMod, only : gridmap_type

  implicit none
  private

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public :: dominant_soil_color
  public :: mkrank

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !

  !EOP
  !===============================================================
contains
  !===============================================================

  !-----------------------------------------------------------------------
  subroutine dominant_soil_color(tgridmap, mask_i, soil_color_i, nsoicol, soil_color_o)
    !
    ! !DESCRIPTION:
    ! Determine the dominant soil color in each output cell
    !
    ! !ARGUMENTS:
    type(gridmap_type) , intent(in)  :: tgridmap
    integer            , intent(in)  :: mask_i(:)       ! input grid: land mask (1 = land, 0 = ocean)
    integer            , intent(in)  :: soil_color_i(:) ! input grid: BATS soil color
    integer            , intent(in)  :: nsoicol         ! number of soil colors
    integer            , intent(out) :: soil_color_o(:) ! output grid: soil color classes
    !
    ! !LOCAL VARIABLES:
    integer, parameter :: num = 2             ! set soil mapunit number
    integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
    integer  :: k, n, ni, no, ns_i, ns_o
    real(r8) :: wt                            ! map overlap weight
    real(r8), allocatable :: wst(:,:)         ! overlap weights, by surface type
    logical :: has_color                      ! whether this grid cell has non-zero color
    integer, parameter :: miss = 99999        ! missing data indicator

    character(len=*), parameter :: subname = 'dominant_soil_color'
    !-----------------------------------------------------------------------

    ns_i = size(mask_i)
    if (size(soil_color_i) /= ns_i) then
       write(6,*) subname, ' ERROR: size of soil_color_i should match size of mask_i'
       write(6,*) 'size(mask_i), size(soil_color_i) = ', &
            size(mask_i), size(soil_color_i)
       call abort()
    end if

    ! find area of overlap for each soil color for each no

    ns_o = size(soil_color_o)
    allocate(wst(0:nsoicol,ns_o))
    wst(0:nsoicol,:) = 0

    ! TODO: need to do a loop to determine
    ! the maximum number of over lap cells throughout the grid 
    ! first get an array that is novr(ns_o) and fill this in - then set
    ! maxovr - to max(novr) - then allocate the array wst to be size of
    ! maxovr,ns_o or 0:nsoilcol,ns_o

    do n = 1,tgridmap%ns
       ni = tgridmap%src_indx(n)
       no = tgridmap%dst_indx(n)
       wt = tgridmap%wovr(n)
       k  = soil_color_i(ni) * mask_i(ni)
       wst(k,no) = wst(k,no) + wt
    enddo

    soil_color_o(:) = 0
    do no = 1,ns_o

       ! If the output cell has any non-zero-colored inputs, then set the weight of
       ! zero-colored inputs to 0, to ensure that the zero-color is NOT dominant.
       if (any(wst(1:nsoicol,no) > 0.)) then
          has_color = .true.
          wst(0,no) = 0.0
       else
          has_color = .false.
       end if

       ! Rank non-zero weights by color type. wsti(1) is the most extensive
       ! color type. 

       if (has_color) then
          call mkrank (nsoicol, wst(0:nsoicol,no), miss, wsti, num)
          soil_color_o(no) = wsti(1)
       end if

       ! If land but no color, set color to 15 (in older dataset generic 
       ! soil color 4)

       if (nsoicol == 8) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 4
          end if
       else if (nsoicol == 20) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 15
          end if
       else
          write(6,*) 'MKSOILCOL error: unhandled nsoicol: ', nsoicol
          call abort()
       end if

       ! Error checks

       if (soil_color_o(no) < 0 .or. soil_color_o(no) > nsoicol) then
          write (6,*) 'MKSOILCOL error: land model soil color = ', &
               soil_color_o(no),' is not valid for lon,lat = ',no
          call abort()
       end if

    end do

    deallocate (wst)

  end subroutine dominant_soil_color


  !-----------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: mkrank
  !
  ! !INTERFACE:
  subroutine mkrank (n, a, miss, iv, num)
    !
    ! !DESCRIPTION:
    ! Return indices of largest [num] values in array [a].
    !
    ! !ARGUMENTS:
    integer , intent(in) :: n        !array length
    real(r8), intent(in) :: a(0:n)   !array to be ranked
    integer , intent(in) :: miss     !missing data value
    integer , intent(in) :: num      !number of largest values requested
    integer , intent(out):: iv(num)  !index to [num] largest values in array [a]
    !
    ! !REVISION HISTORY:
    ! Author: Gordon Bonan
    !
    ! !LOCAL VARIABLES:
    !EOP
    real(r8) a_max       !maximum value in array
    integer i            !array index
    real(r8) delmax      !tolerance for finding if larger value
    integer m            !do loop index
    integer k            !do loop index
    logical exclude      !true if data value has already been chosen
    !-----------------------------------------------------------------------

    delmax = 1.e-06

    ! Find index of largest non-zero number

    iv(1) = miss
    a_max = -9999.

    do i = 0, n
       if (a(i)>0. .and. (a(i)-a_max)>delmax) then
          a_max = a(i)
          iv(1)  = i
       end if
    end do

    ! iv(1) = miss indicates no values > 0. this is an error

    if (iv(1) == miss) then
       write (6,*) 'MKRANK error: iv(1) = missing'
       call abort()
    end if

    ! Find indices of the next [num]-1 largest non-zero number.
    ! iv(m) = miss if there are no more values > 0

    do m = 2, num
       iv(m) = miss
       a_max = -9999.
       do i = 0, n

          ! exclude if data value has already been chosen

          exclude = .false.
          do k = 1, m-1
             if (i == iv(k)) exclude = .true.
          end do

          ! if not already chosen, see if it is the largest of
          ! the remaining values

          if (.not. exclude) then
             if (a(i)>0. .and. (a(i)-a_max)>delmax) then
                a_max = a(i)
                iv(m)  = i
             end if
          end if
       end do
    end do

  end subroutine mkrank

end module mksoilUtilsMod
