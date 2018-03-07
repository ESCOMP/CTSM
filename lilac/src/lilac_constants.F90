module lilac_constants

   implicit none

contains

   integer, parameter :: STRING_128 = 128
   integer, parameter :: STRING_32 = 32
   
   integer, parameter :: FIELD_TYPE_INTEGER = 0
   integer, parameter :: FIELD_TYPE_REAL_8BYTE = 1

   ! known models names
   character(len=*), parameter :: MODEL_NAME_LILAC = 'lilac'
   character(len=*), parameter :: MODEL_NAME_CTSM = 'ctsm'
   character(len=*), parameter :: MODEL_NAME_TEST = 'test'

end module lilac_constants
