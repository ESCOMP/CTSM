CPP = /lib/cpp 
FC  = xlf90_r

.SUFFIXES: .F90 .o

AR = ar
RANLIB	= echo

OBJS = ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
        MeatMod.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o   \
        ESMF_TimeIntervalMod.o ESMF_StubsMod.o ESMF_ModMod.o ESMF_AlarmClockMod.o

libesmf_time.a : $(OBJS)
	$(AR) ru libesmf_time.a $(OBJS)
	$(RANLIB) libesmf_time.a

.F90.o :
	$(CPP) -C -P -I. $*.F90 > $*.f
	$(FC) -c -g $*.f

superclean: 
	/bin/rm -f *.f *.o libesmf_time.a *.mod 

# DEPENDENCIES : only dependencies after this line 

ESMF_AlarmMod.o : ESMF_BaseTimeMod.o ESMF_TimeMod.o ESMF_TimeIntervalMod.o
ESMF_BaseTimeMod.o : ESMF_BaseMod.o
ESMF_ClockMod.o : ESMF_BaseTimeMod.o ESMF_TimeMod.o ESMF_TimeIntervalMod.o
ESMF_AlarmClockMod.o : ESMF_AlarmMod.o ESMF_ClockMod.o
ESMF_TimeMod.o : ESMF_BaseTimeMod.o ESMF_TimeIntervalMod.o
ESMF_BaseMod.o :
ESMF_CalendarMod.o : ESMF_BaseTimeMod.o
ESMF_FractionMod.o : ESMF_BaseTimeMod.o
ESMF_TimeIntervalMod.o : ESMF_BaseTimeMod.o ESMF_CalendarMod.o ESMF_FractionMod.o
ESMF_ModMod.o : ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
        ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o    \
        ESMF_TimeIntervalMod.o MeatMod.o ESMF_StubsMod.o ESMF_AlarmClockMod.o
MeatMod.o : ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
        ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o    \
        ESMF_TimeIntervalMod.o
ESMF_StubsMod.o : ESMF_BaseMod.o ESMF_CalendarMod.o \
        ESMF_AlarmMod.o ESMF_ClockMod.o ESMF_TimeIntervalMod.o ESMF_TimeMod.o

