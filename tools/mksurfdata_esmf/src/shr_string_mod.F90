module shr_string_mod

  implicit none
  private

  public :: shr_string_countChar       ! Count number of char in string, fn
  public :: shr_string_endIndex        ! Index of end of substr in str
  public :: shr_string_betweenTags     ! get the substring between the two tags

!===============================================================================
contains
!===============================================================================

  integer function shr_string_countChar(str,char,rc)

    !  count number of occurances of a single character in a string

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)    ,intent(in)           :: str   ! string to search
    character(1)        ,intent(in)           :: char  ! char to search for
    integer,intent(out),optional :: rc    ! return code

    !----- local -----
    integer :: count    ! counts occurances of char
    integer :: n        ! generic index
    integer :: t01 = 0  ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_countChar) "
    character(*),parameter :: F00     = "('(shr_string_countChar) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    count = 0
    do n = 1, len_trim(str)
       if (str(n:n) == char) count = count + 1
    end do
    shr_string_countChar = count

    if (present(rc)) rc = 0

  end function shr_string_countChar

  !===============================================================================
  subroutine shr_string_betweenTags(string,startTag,endTag,substr,rc)

    !    Get the substring found between the start and end tags.

    ! !INPUT/OUTPUT PARAMETERS:
    character(*)        ,intent(in)  :: string      ! string to search
    character(*)        ,intent(in)  :: startTag    ! start tag
    character(*)        ,intent(in)  :: endTag      ! end tag
    character(*)        ,intent(out) :: substr      ! sub-string between tags
    integer,intent(out),optional :: rc ! retrun code

    !--- local ---
    integer                :: iStart  ! substring start index
    integer                :: iEnd    ! substring end   index
    integer                :: rCode   ! return code
    integer                :: t01 = 0 ! timer
    character(*),parameter :: F00     = "('(shr_string_betweenTags) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    ! * assumes the leading/trailing white space is not part of start & end tags
    !-------------------------------------------------------------------------------

    iStart = shr_string_endIndex(string,trim(adjustL(startTag))) ! end of start tag
    iEnd   =               index(string,trim(adjustL(endTag  ))) ! start of end tag

    rCode = 0
    substr = ""

    if (iStart < 1) then
          write(6,F00) "ERROR: can't find start tag in string"
          write(6,F00) "ERROR: start tag = ",trim(startTag)
          write(6,F00) "ERROR: string    = ",trim(string)
       rCode = 1
    else if (iEnd < 1) then
          write(6,F00) "ERROR: can't find end tag in string"
          write(6,F00) "ERROR: end   tag = ",trim(  endTag)
          write(6,F00) "ERROR: string    = ",trim(string)
       rCode = 2
    else if ( iEnd <= iStart) then
          write(6,F00) "ERROR: start tag not before end tag"
          write(6,F00) "ERROR: start tag = ",trim(startTag)
          write(6,F00) "ERROR: end   tag = ",trim(  endTag)
          write(6,F00) "ERROR: string    = ",trim(string)
       rCode = 3
    else if ( iStart+1 == iEnd ) then
       substr = ""
       write(6,F00) "WARNING: zero-length substring found in ",trim(string)
    else
       substr = string(iStart+1:iEnd-1)
       if (len_trim(substr) == 0 ) then
          write(6,F00) "WARNING: white-space substring found in ",trim(string)
       end if
    end if

    if (present(rc)) rc = rCode

  end subroutine shr_string_betweenTags

  !===============================================================================
  integer function shr_string_endIndex(string,substr,rc)

    !  Get the ending index of substr within string

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) ,intent(in)           :: string ! string to search
    character(len=*) ,intent(in)           :: substr ! sub-string to search for
    integer          ,intent(out),optional :: rc     ! return code

    !--- local ---
    integer   :: i       ! generic index
    !-------------------------------------------------------------------------------
    ! Notes:
    ! * returns zero if substring not found, uses len_trim() intrinsic
    ! * very similar to: i = index(str,substr,back=.true.) 
    ! * do we need this function?
    !-------------------------------------------------------------------------------

    i = index(trim(string),trim(substr))
    if ( i == 0 ) then
       shr_string_endIndex = 0  ! substr is not in string
    else
       shr_string_endIndex = i + len_trim(substr) - 1
    end if

    !  -------------------------------------------------------------------
    !  i = index(trim(string),trim(substr),back=.true.)
    !  if (i == len(string)+1) i = 0
    !  shr_string_endIndex = i
    !  -------------------------------------------------------------------

    if (present(rc)) rc = 0

  end function shr_string_endIndex

end module shr_string_mod
