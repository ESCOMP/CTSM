# Pick up the shared perf_timers_mod from the parent perf_testing/ dir.
VPATH := ..

OBJ := SoilBiogeochemCompetition.o perf_timers_mod.o driver.o

# Module-use ordering. driver.F90 uses SoilBiogeochemCompetition_mod;
# both driver.F90 and SoilBiogeochemCompetition.F90 use perf_timers_mod
# (start/stop calls inside the routine; print/dump calls in the driver).
driver.o: SoilBiogeochemCompetition.o perf_timers_mod.o
SoilBiogeochemCompetition.o: perf_timers_mod.o

include ../Makefile.common
