FC      ?= nvfortran
FFLAGS  ?= -O3 -g
LDFLAGS ?=

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
