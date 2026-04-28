FC      ?= nvfortran
FFLAGS  ?= -O3 -g
LDFLAGS ?=

# TIMING=1 (default) compiles the driver with built-in system_clock-based
# timing. TIMING=0 omits it, leaving wall-clock measurement to an external
# profiler.
TIMING  ?= 1
ifeq ($(TIMING),1)
    FFLAGS += -DPERF_TIMING
endif

OBJ := p2c_2d_filter.o driver.o

driver: $(OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ)

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

# driver.F90 uses the p2c_2d_filter_mod module
driver.o: p2c_2d_filter.o

.PHONY: clean
clean:
	./clean.sh
