      integer GPTLwall
      integer GPTLcpu
      integer GPTLabort_on_error
      integer GPTLoverhead
      integer GPTLdepthlimit

      parameter (GPTLwall           = 1)
      parameter (GPTLcpu            = 2)
      parameter (GPTLabort_on_error = 3)
      parameter (GPTLoverhead       = 4)
      parameter (GPTLdepthlimit     = 5)

! Externals: new interface.  setoption and initialize are now functions

      integer gptlsetoption, gptlinitialize
      external gptlsetoption
      external gptlinitialize
!      external gptlstart
!      external gptlstop
!      external gptlstamp 
!      external gptlpr
!      external gptlreset 
!      external gptlfinalize

! Externals: old interface.  All entries were subroutines

      external t_setoptionf
      external t_initializef
      external t_startf
      external t_stopf
      external t_stampf 
      external t_prf
      external t_resetf 
      external t_finalizef
