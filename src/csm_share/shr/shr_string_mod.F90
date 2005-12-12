!===============================================================================
! CVS: $Id$
! CVS: $Source$
! CVS: $Name$
!===============================================================================
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_string_mod -- string and list methods
!
! !DESCRIPTION:
!    General string and specific list method.  A list is a single string
!    that is delimited by a character forming multiple fields, ie,
!    character(len=*) :: mylist = "t:s:u1:v1:u2:v2:taux:tauy"
!    The delimiter is called listDel in this module, is default ":",
!    but can be set by a call to shr_string_listSetDel.
!
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

module shr_string_mod

! !USES:

   use shr_kind_mod   ! F90 kinds
   use shr_sys_mod    ! shared system calls
   use shr_cal_mod    ! shared calendar

   implicit none
   private

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_string_countChar       ! Count number of char in string, fn
   public :: shr_string_lastIndex       ! Index of last substr in str
   public :: shr_string_endIndex        ! Index of end of substr in str
   public :: shr_string_leftAlign       ! remove leading white space
   public :: shr_string_alphanum        ! remove all non alpha-numeric characters
   public :: shr_string_betweenTags     ! get the substring between the two tags
   public :: shr_string_parseCFtunit    ! parse CF time units
   public :: shr_string_clean           ! Set string to all white space
   public :: shr_string_listIsValid     ! test for a valid "list"
   public :: shr_string_listGetNum      ! Get number of fields in list, fn
   public :: shr_string_listGetIndex    ! Get index of field
   public :: shr_string_listGetIndexF   ! function version of listGetIndex
   public :: shr_string_listGetName     ! get k-th field name
   public :: shr_string_listIntersect   ! get intersection of two field lists
   public :: shr_string_listUnion       ! get union of two field lists
   public :: shr_string_listMerge       ! merge two lists to form third
   public :: shr_string_listAppend      ! append list at end of another
   public :: shr_string_listPrepend     ! prepend list in front of another
   public :: shr_string_listSetDel      ! Set field delimeter in lists
   public :: shr_string_listGetDel      ! Get field delimeter in lists
   public :: shr_string_setAbort        ! set local abort flag
   public :: shr_string_setDebug        ! set local debug flag

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP

   character(len=1)    ,save :: listDel = ":"    ! note single exec implications
   logical             ,save :: doabort = .true.
   integer(SHR_KIND_IN),save :: debug   = 0

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_countChar -- Count number of occurances of a character
!
! !DESCRIPTION:
!  count number of occurances of a single character in a string
!     \newline
!     n = shr\_string\_countChar(string,character)
!
! !REVISION HISTORY:
!     2005-Feb-28 - First version from dshr_bundle
!
! !INTERFACE: ------------------------------------------------------------------

integer function shr_string_countChar(str,char,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: str   ! string to search
   character(1)        ,intent(in)           :: char  ! char to search for
   integer(SHR_KIND_IN),intent(out),optional :: rc    ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: count    ! counts occurances of char
   integer(SHR_KIND_IN) :: n        ! generic index

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
!BOP ===========================================================================
!
! !IROUTINE: shr_string_lastIndex -- Get index of last substr within string
!
! !DESCRIPTION:
!  Get index of last substr within string
!     \newline
!     n = shr\_string\_lastIndex(string,substring)
!
! !REVISION HISTORY:
!     2005-Feb-28 - First version from dshr_domain
!
! !INTERFACE: ------------------------------------------------------------------

integer function shr_string_lastIndex(string,substr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: string ! string to search
   character(*)        ,intent(in)           :: substr ! sub-string to search for
   integer(SHR_KIND_IN),intent(out),optional :: rc     ! return code

!EOP

   !--- local ---
   integer(SHR_KIND_IN)   :: i       ! generic index
   character(SHR_KIND_CL) :: str     ! temporary work string
   character(1)           :: char    ! temporary work string

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_lastIndex) "
   character(*),parameter :: F00     = "('(shr_string_lastIndex) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   str=' '
   str=string

   if (substr(1:1) == 'x') char='y'
   if (substr(1:1) /= 'x') char='x'

   shr_string_lastIndex = index(str,substr)
   i = shr_string_lastIndex
   do while (i > 0)
      !--- search for another occurance of substr ---
      shr_string_lastIndex = i
      str(i:i) = char ! so we won't find this instance of substr again
      i = index(str,substr)
   end do

   if (present(rc)) rc = 0

end function shr_string_lastIndex

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_endIndex -- Get the ending index of substr within string
!
! !DESCRIPTION:
!  Get the ending index of substr within string
!     \newline
!     n = shr\_string\_endIndex(string,substring)
!
! !REVISION HISTORY:
!     2005-May-10 - B. Kauffman, first version.
!
! !INTERFACE: ------------------------------------------------------------------

integer function shr_string_endIndex(string,substr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: string ! string to search
   character(*)        ,intent(in)           :: substr ! sub-string to search for
   integer(SHR_KIND_IN),intent(out),optional :: rc     ! return code

!EOP

   !--- local ---
   integer(SHR_KIND_IN)   :: i       ! generic index

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_endIndex) "
   character(*),parameter :: F00     = "('(shr_string_endIndex) ',4a)"

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

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_leftAlign -- remove leading white space
!
! !DESCRIPTION:
!    Remove leading white space
!     \newline
!     call shr\_string\_leftAlign(string)
!
! !REVISION HISTORY:
!     2005-Apr-28 - B. Kauffman - First version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_leftAlign(str,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(inout)          :: str
   integer(SHR_KIND_IN),intent(out)  ,optional :: rc   ! return code

!EOP

   !----- local ----
   integer(SHR_KIND_IN) :: rCode ! return code

!-------------------------------------------------------------------------------
! note: 
! * ?? this routine isn't needed, use the intrisic adjustL instead ??
!-------------------------------------------------------------------------------

!  -------------------------------------------------------------------
!  --- I used this until I discovered the intrinsic function below - BK
!  do while (len_trim(str) > 0 ) 
!     if (str(1:1) /= ' ') exit
!     str = str(2:len_trim(str))
!  end do
!  rCode = 0
!  !! (len_trim(str) == 0 ) rCode = 1  ! ?? appropriate ??
!  -------------------------------------------------------------------

   str = adjustL(str)
   if (present(rc)) rc = 0

end subroutine shr_string_leftAlign

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_alphanum -- remove non alpha numeric characters
!
! !DESCRIPTION:
!    Remove all non alpha numeric characters from string
!     \newline
!     call shr\_string\_alphanum(string)
!
! !REVISION HISTORY:
!     2005-Aug-01 - T. Craig - First version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_alphanum(str,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(inout)          :: str
   integer(SHR_KIND_IN),intent(out)  ,optional :: rc   ! return code

!EOP

   !----- local ----
   integer(SHR_KIND_IN) :: rCode  ! return code
   integer(SHR_KIND_IN) :: n,icnt ! counters

!-------------------------------------------------------------------------------

   icnt = 0
   do n=1,len_trim(str)
     if ((str(n:n) >= 'a' .and. str(n:n) <= 'z') .or.  &
         (str(n:n) >= 'A' .and. str(n:n) <= 'Z') .or.  &
         (str(n:n) >= '0' .and. str(n:n) <= '9')) then
       icnt = icnt + 1
       str(icnt:icnt) = str(n:n)
     endif
   enddo
   do n=icnt+1,len(str)
     str(n:n) = ' '
   enddo

   if (present(rc)) rc = 0

end subroutine shr_string_alphanum

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_betweenTags -- Get the substring between the two tags.
!
! !DESCRIPTION:
!    Get the substring found between the start and end tags.
!    \newline
!    call shr\_string\_betweenTags(string,startTag,endTag,substring,rc)
!
! !REVISION HISTORY:
!     2005-May-11 - B. Kauffman, first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_betweenTags(string,startTag,endTag,substr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)  :: string      ! string to search
   character(*)        ,intent(in)  :: startTag    ! start tag
   character(*)        ,intent(in)  :: endTag      ! end tag
   character(*)        ,intent(out) :: substr      ! sub-string between tags
   integer(SHR_KIND_IN),intent(out),optional :: rc ! retrun code

!EOP

   !--- local ---
   integer(SHR_KIND_IN)   :: iStart  ! substring start index
   integer(SHR_KIND_IN)   :: iEnd    ! substring end   index
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_betweenTags) "
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
      if (len_trim(substr) == 0 ) &
         & write(6,F00) "WARNING: white-space substring found in ",trim(string)
   end if

   if (present(rc)) rc = rCode

end subroutine shr_string_betweenTags

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_parseCFtunit -- Parse CF time unit
!
! !DESCRIPTION:
!  Parse CF time unit into a delta string name and a base time in yyyymmdd
!  and seconds (nearest integer actually).
!     \newline
!     call shr\_string\_parseCFtunit(string,substring)
!     \newline
!  Input string is like "days since 0001-06-15 15:20:45.5 -6:00"
!    - recognizes "days", "hours", "minutes", "seconds"
!    - must have at least yyyy-mm-dd, hh:mm:ss.s is optional
!    - expects a "since" in the string
!    - ignores time zone part
!
! !REVISION HISTORY:
!     2005-May-15 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_parseCFtunit(string,unit,bdate,bsec,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: string ! string to search
   character(*)        ,intent(out)          :: unit   ! delta time unit
   integer(SHR_KIND_IN),intent(out)          :: bdate  ! base date yyyymmdd
   real(SHR_KIND_R8)   ,intent(out)          :: bsec   ! base seconds
   integer(SHR_KIND_IN),intent(out),optional :: rc     ! return code

!EOP

   !--- local ---
   integer(SHR_KIND_IN)   :: i,i1,i2 ! generic index
   character(SHR_KIND_CL) :: tbase   ! baseline time
   character(SHR_KIND_CL) :: lstr    ! local string
   integer(SHR_KIND_IN)   :: yr,mo,da,hr,min  ! time stuff
   real(SHR_KIND_R8)      :: sec              ! time stuff

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_parseCFtunit) "
   character(*),parameter :: F00     = "('(shr_string_parseCFtunit) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   unit = 'none'
   bdate = 0
   bsec = 0.0_SHR_KIND_R8

   i = shr_string_lastIndex(string,'days ')
   if (i > 0) unit = 'days'
   i = shr_string_lastIndex(string,'hours ')
   if (i > 0) unit = 'hours'
   i = shr_string_lastIndex(string,'minutes ')
   if (i > 0) unit = 'minutes'
   i = shr_string_lastIndex(string,'seconds ')
   if (i > 0) unit = 'seconds'

   if (trim(unit) == 'none') then
     write(6,F00) ' ERROR time unit unknown'
     call shr_string_abort(subName//' time unit unknown')
   endif

   i = shr_string_lastIndex(string,' since ')
   if (i < 1) then
     write(6,F00) ' ERROR since does not appear in unit attribute for time '
     call shr_string_abort(subName//' no since in attr name')
   endif
   tbase = trim(string(i+6:))
   call shr_string_leftAlign(tbase)

   if (debug > 0) then
     write(6,*) trim(subName)//' '//'unit '//trim(unit)
     write(6,*) trim(subName)//' '//'tbase '//trim(tbase)
   endif

   yr=0; mo=0; da=0; hr=0; min=0; sec=0
   i1 = 1

   i2 = index(tbase,'-') - 1
   lstr = tbase(i1:i2)
   read(lstr,*,ERR=200,END=200) yr
   tbase = tbase(i2+2:)
   call shr_string_leftAlign(tbase)

   i2 = index(tbase,'-') - 1
   lstr = tbase(i1:i2)
   read(lstr,*,ERR=200,END=200) mo
   tbase = tbase(i2+2:)
   call shr_string_leftAlign(tbase)

   i2 = index(tbase,' ') - 1
   lstr = tbase(i1:i2)
   read(lstr,*,ERR=200,END=200) da
   tbase = tbase(i2+2:)
   call shr_string_leftAlign(tbase)

   i2 = index(tbase,':') - 1
   lstr = tbase(i1:i2)
   read(lstr,*,ERR=200,END=100) hr
   tbase = tbase(i2+2:)
   call shr_string_leftAlign(tbase)

   i2 = index(tbase,':') - 1
   lstr = tbase(i1:i2)
   read(lstr,*,ERR=200,END=100) min
   tbase = tbase(i2+2:)
   call shr_string_leftAlign(tbase)

   i2 = index(tbase,' ') - 1
   lstr = tbase(i1:i2)
   read(lstr,*,ERR=200,END=100) sec

100  continue

   if (debug > 0) write(6,*) trim(subName),'ymdhms:',yr,mo,da,hr,min,sec

   call shr_cal_ymd2date(yr,mo,da,bdate)
   bsec = real(hr*3600 + min*60,SHR_KIND_R8) + sec

   if (present(rc)) rc = 0

   return

200  write(6,F00) 'ERROR 200 on char num read '
     call shr_string_abort(subName//' ERROR on char num read')
     return

end subroutine shr_string_parseCFtunit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_clean -- Clean a string, set it to "blank"
!
! !DESCRIPTION:
!     Clean a string, set it to blank
!     \newline
!     call shr\_string\_clean(string,rc)
!
! !REVISION HISTORY:
!     2005-May-05 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_clean(string,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(inout) :: string  ! list/string
   integer(SHR_KIND_IN),optional,intent(out)   :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: n       ! counter
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_clean) "
   character(*),parameter :: F00     = "('(shr_string_clean) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0
   do n=1,len(string)
     string(n:n) = ' '
   enddo
   if (present(rc)) rc = rCode

end subroutine shr_string_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listIsValid -- determine whether string is a valid list
!
! !DESCRIPTION:
!     Determine whether string is a valid list
!     \newline
!     logical_var = shr\_string\_listIsValid(list,rc)
!
! !REVISION HISTORY:
!     2005-May-05 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

logical function shr_string_listIsValid(list,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(in)  :: list    ! list/string
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer  (SHR_KIND_IN) :: i,k     ! generic loop index
   integer  (SHR_KIND_IN) :: kFlds   ! number of fields in list
   character(SHR_KIND_CL) :: subList ! sub-string/list
   integer  (SHR_KIND_IN) :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listIsValid) "
   character(*),parameter :: F00     = "('(shr_string_listIsValid) ',4a)"

!-------------------------------------------------------------------------------
! check that each field in list has at least one character
!-------------------------------------------------------------------------------

   rCode = 0

   kFlds = shr_string_listGetNum(list)
   subList = list
   do k=1,kFlds-1
      i = index(subList,listDel)               ! position  of 1st delimiter
      if (i < 2) rCode = 1
      subList = subList(i+1:len_trim(subList)) ! delete up to 1st delimiter
   end do
   if (len_trim(subList) == 0) rCode = 1       ! n=kFlds ~ last field in list
   
   if (rCode == 0) shr_string_listIsValid = .true.
   if (rCode /= 0) shr_string_listIsValid = .false.
   if (present(rc)) rc = rCode

end function shr_string_listIsValid

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listGetName -- Get name of k-th field in list
!
! !DESCRIPTION:
!     Get name of k-th field in list
!     \newline
!     call shr\_string\_listGetName(list,k,name,rc)
!
! !REVISION HISTORY:
!     2005-May-05 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listGetName(list,k,name,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(in)  :: list    ! list/string
   integer(SHR_KIND_IN)         ,intent(in)  :: k       ! index of field
   character(*)                 ,intent(out) :: name    ! k-th name in list
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: i,j,n   ! generic indecies
   integer(SHR_KIND_IN)   :: kFlds   ! number of fields in list
   character(SHR_KIND_CL) :: str     ! temp string
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listGetName) "
   character(*),parameter :: F00     = "('(shr_string_listGetName) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   !--- check that this is a valid list ---
   if (.not. shr_string_listIsValid(list,rCode) ) then
      write(6,F00) "ERROR: invalid list = ",trim(list)
      call shr_string_abort(subName//" ERROR: invalid list = "//trim(list))
   end if

   !--- check that this is a valid index ---
   kFlds = shr_string_listGetNum(list)
   if (k<1 .or. kFlds<k) then
      write(6,*) subName,"ERROR: invalid index = ",k
      write(6,*) subName,"ERROR:          list = ",trim(list)
      call shr_string_abort(subName//" ERROR: invalid index")
   end if

   !--- start with whole list, then remove fields before and after desired field ---
   str = list

   !--- remove field names before desired field ---
   do n=2,k
      i = index(str,listDel)
      str = str(i+1:len_trim(str))
   end do

   !--- remove field names after desired field ---
   if ( k < kFlds ) then
      j = index(str,listDel)
      str = str(1:j-1)
   end if

   !--- copy result into output variable ---
   name = str

   if (present(rc)) rc = rCode

end subroutine shr_string_listGetName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listIntersect -- Get intersection of two field lists
!
! !DESCRIPTION:
!     Get intersection of two fields lists, write into third list
!     \newline
!     call shr\_string\_listIntersect(list1,list2,listout)
!
! !REVISION HISTORY:
!     2005-May-05 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listIntersect(list1,list2,listout,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(in)  :: list1   ! list/string
   character(*)                 ,intent(in)  :: list2   ! list/string
   character(*)                 ,intent(out) :: listout ! list/string
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: nf,n1,n2 ! counters
   character(SHR_KIND_CS) :: name     ! field name
   integer(SHR_KIND_IN)   :: rCode    ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listIntersect) "
   character(*),parameter :: F00     = "('(shr_string_listIntersect) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   nf = shr_string_listGetNum(list1)
   call shr_string_clean(listout)
   do n1 = 1,nf
     call shr_string_listGetName(list1,n1,name,rCode)
     n2 = shr_string_listGetIndexF(list2,name)
     if (n2 > 0) then
       call shr_string_listAppend(listout,name)
     endif
   enddo

   if (present(rc)) rc = rCode

end subroutine shr_string_listIntersect

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listUnion -- Get union of two field lists
!
! !DESCRIPTION:
!     Get union of two fields lists, write into third list
!     \newline
!     call shr\_string\_listUnion(list1,list2,listout)
!
! !REVISION HISTORY:
!     2005-May-05 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listUnion(list1,list2,listout,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(in)  :: list1   ! list/string
   character(*)                 ,intent(in)  :: list2   ! list/string
   character(*)                 ,intent(out) :: listout ! list/string
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: nf,n1,n2 ! counters
   character(SHR_KIND_CS) :: name     ! field name
   integer(SHR_KIND_IN)   :: rCode    ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listUnion) "
   character(*),parameter :: F00     = "('(shr_string_listUnion) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   call shr_string_clean(listout)

   nf = shr_string_listGetNum(list1)
   do n1 = 1,nf
     call shr_string_listGetName(list1,n1,name,rCode)
     n2 = shr_string_listGetIndexF(listout,name)
     if (n2 < 1) then
       call shr_string_listAppend(listout,name)
     endif
   enddo

   nf = shr_string_listGetNum(list2)
   do n1 = 1,nf
     call shr_string_listGetName(list2,n1,name,rCode)
     n2 = shr_string_listGetIndexF(listout,name)
     if (n2 < 1) then
       call shr_string_listAppend(listout,name)
     endif
   enddo

   if (present(rc)) rc = rCode

end subroutine shr_string_listUnion

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listMerge -- Merge lists two list to third
!
! !DESCRIPTION:
!     Merge two list to third
!     \newline
!     call shr\_string\_listMerge(list1,list2,listout)
!     call shr\_string\_listMerge(list1,list2,list1)
!
! !REVISION HISTORY:
!     2005-May-05 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listMerge(list1,list2,listout,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(in)  :: list1   ! list/string
   character(*)                 ,intent(in)  :: list2   ! list/string
   character(*)                 ,intent(out) :: listout ! list/string
   integer(SHR_KIND_IN),optional,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   character(SHR_KIND_CL) :: l1,l2   ! local char strings
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listMerge) "
   character(*),parameter :: F00     = "('(shr_string_listMerge) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   call shr_string_clean(l1)
   call shr_string_clean(l2)
   call shr_string_clean(listout)
   if (len_trim(list1) > len(l1)) call shr_string_abort(subName//' ERROR l1 len')
   if (len_trim(list2) > len(l2)) call shr_string_abort(subName//' ERROR l2 len')
   l1 = trim(list1)
   l2 = trim(list2)
   call shr_string_leftAlign(l1,rCode)
   call shr_string_leftAlign(l2,rCode)
   if (len_trim(l1)+len_trim(l2)+1 > len(listout)) &
     call shr_string_abort(subName//' ERROR listout len')
   if (len_trim(l1) == 0) then
     listout = trim(l2)
   else
     listout = trim(l1)//":"//trim(l2)
   endif

   if (present(rc)) rc = rCode

end subroutine shr_string_listMerge

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listAppend -- Append one list to another
!
! !DESCRIPTION:
!     Append one list to another
!     \newline
!     call shr\_string\_listAppend(list,listadd)
!
! !REVISION HISTORY:
!     2005-May-05 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listAppend(list,listadd,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(inout) :: list   ! list/string
   character(*)                 ,intent(in)    :: listadd ! list/string
   integer(SHR_KIND_IN),optional,intent(out)   :: rc      ! return code

!EOP

   !----- local -----
   character(SHR_KIND_CL) :: l1      ! local string
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listAppend) "
   character(*),parameter :: F00     = "('(shr_string_listAppend) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   if (len_trim(listadd) > len(l1)) call shr_string_abort(subName//' ERROR l1 len')
   call shr_string_clean(l1)
   l1 = trim(listadd)
   call shr_string_leftAlign(l1,rCode)
   if (len_trim(list)+len_trim(l1)+1 > len(list)) &
     call shr_string_abort(subName//' ERROR list len')
   if (len_trim(list) == 0) then
     list = trim(l1)
   else
     list = trim(list)//":"//trim(l1)
   endif

   if (present(rc)) rc = rCode

end subroutine shr_string_listAppend

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listPrepend -- Prepend one list to another
!
! !DESCRIPTION:
!     Prepend one list to another
!     \newline
!     call shr\_string\_listPrepend(listadd,list)
!     \newline
!     results in listadd:list
!
! !REVISION HISTORY:
!     2005-May-05 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listPrepend(listadd,list,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                 ,intent(in)    :: listadd ! list/string
   character(*)                 ,intent(inout) :: list   ! list/string
   integer(SHR_KIND_IN),optional,intent(out)   :: rc      ! return code

!EOP

   !----- local -----
   character(SHR_KIND_CL) :: l1      ! local string
   integer(SHR_KIND_IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listPrepend) "
   character(*),parameter :: F00     = "('(shr_string_listPrepend) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   rCode = 0

   if (len_trim(listadd) > len(l1)) call shr_string_abort(subName//' ERROR l1 len')
   call shr_string_clean(l1)
   l1 = trim(listadd)
   call shr_string_leftAlign(l1,rCode)
   call shr_string_leftAlign(list,rCode)
   if (len_trim(list)+len_trim(l1)+1 > len(list)) &
     call shr_string_abort(subName//' ERROR list len')
   if (len_trim(l1) == 0) then
     list = trim(list)
   else
     list = trim(l1)//":"//trim(list)
   endif

   if (present(rc)) rc = rCode

end subroutine shr_string_listPrepend

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listGetIndexF -- Get index of field in string
!
! !DESCRIPTION:
!     Get index of field in string
!     \newline
!     k = shr\_string\_listGetIndex(str,"taux")
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

integer function shr_string_listGetIndexF(string,fldStr)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in) :: string   ! string
   character(*),intent(in) :: fldStr   ! name of field

!EOP

   !----- local -----
   integer(SHR_KIND_IN)    :: k        ! local index variable
   integer(SHR_KIND_IN)    :: rc       ! error code

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listGetIndexF) "
   character(*),parameter :: F00     = "('(shr_string_listGetIndexF) ',4a)"

!-------------------------------------------------------------------------------

   call shr_string_listGetIndex(string,fldStr,k,print=.false.,rc=rc)
   shr_string_listGetIndexF = k

end function shr_string_listGetIndexF

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listGetIndex -- Get index of field in string
!
! !DESCRIPTION:
!     Get index of field in string
!     \newline
!     call shr\_string\_listGetIndex(str,"taux",k,rc)
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listGetIndex(string,fldStr,k,print,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: string  ! string
   character(*)        ,intent(in)           :: fldStr  ! name of field
   integer(SHR_KIND_IN),intent(out)          :: k       ! index of field
   logical             ,intent(in) ,optional :: print   ! print switch
   integer(SHR_KIND_IN),intent(out),optional :: rc      ! return code

!EOP

   !----- local -----
   character(SHR_KIND_CL) :: fieldNames       ! string of field names
   character(SHR_KIND_CL) :: str              ! string with one field name
   integer(SHR_KIND_IN)   :: n                ! index for colon position
   integer(SHR_KIND_IN)   :: nFields          ! number of fields in a string
   logical                :: found            ! T => field found in fieldNames
   logical                :: lprint           ! local print flag

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listGetIndex) "
   character(*),parameter :: F00     = "('(shr_string_listGetIndex) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   lprint = .false.
   if (present(print)) lprint = print

   !--- confirm proper size of input data ---
   if (len_trim(trim(fldStr)) < 1) then
      if (lprint) write(6,F00) "ERROR: input field name has 0 length"
      call shr_string_abort(subName//"invalid field name")
   end if

   !--- search for field name in string's list of fields ---
   found      = .false.
   fieldNames = string
   nFields    = shr_string_listGetNum(string)
   do k = 1, nFields
      n = index(fieldNames,listDel)
      !--- sanity check ---
      if ((k <nFields .and. n<1) .or. (k==nFields .and. n>0)) then
	 call shr_string_abort(subName//"ERROR: wrong string%nf ?")
      end if
      !--- get field name of field number k ---
      if ( n == 0) then   
         str = fieldNames        ! this is the last field name in fieldNames
      else               
         str = fieldNames(1:n-1) ! *not* the last field name in fieldNames
      endif
      !--- is it a match? ---
      if (trim(fldStr) == trim(str)) found = .true.
      if (found) exit
      !--- remove 1st field & repeat ---
      fieldNames = fieldNames(n+1:len_trim(fieldNames)) 
   end do

   !--- not finding a field is not a fatal error ---
   if (.not. found) then
      k = 0
      if (lprint) write(6,F00) "FYI: field ",trim(fldStr)," not found in list ",trim(string)
      if (present(rc)) rc = 1
   end if

end subroutine shr_string_listGetIndex

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listGetNum -- get number of fields in a string list
!
! !DESCRIPTION:
!  return number of fields in string list
!
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - First version
!
! !INTERFACE: ------------------------------------------------------------------

integer function shr_string_listGetNum(str)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in) :: str   ! string to search

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: count    ! counts occurances of char

   !----- formats -----
   character(*),parameter :: subName =   "(shr_string_listGetNum) "
   character(*),parameter :: F00     = "('(shr_string_listGetNum) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   count = shr_string_countChar(str,listDel)
   shr_string_listGetNum = count + 1

end function shr_string_listGetNum

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listSetDel -- Set list delimeter character
!
! !DESCRIPTION:
!     Set field delimeter character in lists
!     \newline
!     call shr\_string\_listSetDel(":")
!
! !REVISION HISTORY:
!     2005-Apr-30  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listSetDel(cflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(len=1),intent(in) :: cflag

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_string_listSetDel) "
  character(*),parameter :: F00     = "('(shr_string_listSetDel) ',a) "

!-------------------------------------------------------------------------------

  if (debug > 0) write(6,F00) 'changing listDel from '//trim(listDel)//' to '//trim(cflag)
  listDel = trim(cflag)

end subroutine shr_string_listSetDel

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_listGetDel -- Get list delimeter character
!
! !DESCRIPTION:
!     Get field delimeter character in lists
!     \newline
!     call shr\_string\_listGetDel(del)
!
! !REVISION HISTORY:
!     2005-May-15  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_listGetDel(del)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),intent(out) :: del

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_string_listGetDel) "
  character(*),parameter :: F00     = "('(shr_string_listGetDel) ',a) "

!-------------------------------------------------------------------------------

  del = trim(listDel)

end subroutine shr_string_listGetDel

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_setAbort -- Set local shr_string abort flag
!
! !DESCRIPTION:
!     Set local shr_string abort flag, true = abort, false = print and continue
!     \newline
!     call shr\_string\_setAbort(.false.)
!
! !REVISION HISTORY:
!     2005-Apr-30  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical,intent(in) :: flag

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_string_setAbort) "
  character(*),parameter :: F00     = "('(shr_string_setAbort) ',a) "

!-------------------------------------------------------------------------------

  if (debug > 0) then
    if (flag) then
      write(6,F00) 'setting debug to true'
    else
      write(6,F00) 'setting debug to false'
    endif
  endif

  doabort = flag

end subroutine shr_string_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_string_setDebug -- Set local shr_string debug level
!
! !DESCRIPTION:
!     Set local shr_string debug level, 0 = production
!     \newline
!     call shr\_string\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-Apr-30  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_string_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_string_setDebug) "
  character(*),parameter :: F00     = "('(shr_string_setDebug) ',a) "
  character(*),parameter :: F01     = "('(shr_string_setDebug) ',a,i3,a,i3) "

!-------------------------------------------------------------------------------

  if (debug > 0.or.iflag > 0) write(6,F01) 'changing debug level from ',debug,' to ',iflag
  debug = iflag

end subroutine shr_string_setDebug

!===============================================================================
!===============================================================================

subroutine shr_string_abort(string)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(in) :: string

!EOP

  !--- local ---
  character(SHR_KIND_CL) :: lstring
  character(*),parameter :: subName =   "(shr_string_abort)"
  character(*),parameter :: F00     = "('(shr_string_abort) ',a)"

!-------------------------------------------------------------------------------

  lstring = ''
  if (present(string)) lstring = string

  if (doabort) then
    call shr_sys_abort(lstring)
  else
    write(6,F00) ' no abort:'//trim(lstring)
  endif

end subroutine shr_string_abort

!===============================================================================
!===============================================================================

end module shr_string_mod
