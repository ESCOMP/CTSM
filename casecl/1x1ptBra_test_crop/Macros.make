PIO_FILESYSTEM_HINTS := gpfs
PNETCDF_PATH := $(PNETCDF)
NETCDF_PATH := $(NETCDF)
SUPPORTS_CXX := FALSE
ifeq ($(COMPILER),intel)
  MPIFC :=  mpif90 
  FFLAGS_NOOPT :=  -O0 
  MPICC :=  mpicc  
  SCC :=  icc 
  MPICXX :=  mpicxx 
  HAS_F2008_CONTIGUOUS := FALSE
  CXX_LDFLAGS :=  -cxxlib 
  SUPPORTS_CXX := TRUE
  FFLAGS :=  -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source  
  FIXEDFLAGS :=  -fixed -132 
  CXX_LINKER := FORTRAN
  FC_AUTO_R8 :=  -r8 
  CFLAGS :=   -qno-opt-dynamic-align -fp-model precise -std=gnu99 
  FREEFLAGS :=  -free 
  SFC :=  ifort 
  SCXX :=  icpc 
  ifeq ($(MPILIB),mpt)
    ifeq ($(compile_threaded),true)
      PFUNIT_PATH := $(CESMDATAROOT)/tools/pFUnit/pFUnit3.2.8_cheyenne_Intel17.0.1_MPI_openMP
    endif
  endif
  ifeq ($(MPILIB),mpi-serial)
    ifeq ($(compile_threaded),false)
      PFUNIT_PATH := $(CESMDATAROOT)/tools/pFUnit/pFUnit3.2.8_cheyenne_Intel17.0.1_noMPI_noOpenMP
    endif
  endif
endif
ifeq ($(COMPILER),gnu)
  MPIFC :=  mpif90 
  FFLAGS_NOOPT :=  -O0 
  MPICC :=  mpicc  
  SCC :=  gcc 
  MPICXX :=  mpicxx 
  SUPPORTS_CXX := TRUE
  FFLAGS :=   -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none 
  FIXEDFLAGS :=   -ffixed-form 
  CXX_LINKER := FORTRAN
  FC_AUTO_R8 :=  -fdefault-real-8 
  CFLAGS :=  -std=gnu99 
  FREEFLAGS :=  -ffree-form 
  SFC :=  gfortran 
  SCXX :=  g++ 
endif
CPPDEFS := $(CPPDEFS)  -D$(OS) 
ifeq ($(MODEL),gptl)
  CPPDEFS := $(CPPDEFS)  -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY 
endif
ifeq ($(MODEL),pop)
  CPPDEFS := $(CPPDEFS)  -D_USE_FLOW_CONTROL 
endif
ifeq ($(COMPILER),intel)
  FFLAGS := $(FFLAGS)  -qopt-report -xCORE_AVX2 -no-fma
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DCPRINTEL
  CFLAGS := $(CFLAGS)  -qopt-report -xCORE_AVX2 -no-fma
  ifeq ($(compile_threaded),true)
    FFLAGS := $(FFLAGS)  -qopenmp 
    CFLAGS := $(CFLAGS)  -qopenmp 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created 
    CMAKE_OPTS := $(CMAKE_OPTS)  -DPIO_ENABLE_LOGGING=ON 
    CFLAGS := $(CFLAGS)  -O0 -g 
  endif
  ifeq ($(DEBUG),FALSE)
    FFLAGS := $(FFLAGS)  -O2 -debug minimal 
    CFLAGS := $(CFLAGS)  -O2 -debug minimal 
  endif
  ifeq ($(MPILIB),mpt)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),openmpi)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpi-serial)
    SLIBS := $(SLIBS)  -mkl 
  endif
  ifeq ($(compile_threaded),true)
    FFLAGS_NOOPT := $(FFLAGS_NOOPT)  -qopenmp 
    LDFLAGS := $(LDFLAGS)  -qopenmp 
  endif
endif
ifeq ($(COMPILER),gnu)
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU
  ifeq ($(compile_threaded),true)
    FFLAGS := $(FFLAGS)  -fopenmp 
    CFLAGS := $(CFLAGS)  -fopenmp 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -g -Wall -Og -fbacktrace -ffpe-trap=zero,overflow -fcheck=bounds 
    CFLAGS := $(CFLAGS)  -g -Wall -Og -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds 
  endif
  ifeq ($(DEBUG),FALSE)
    FFLAGS := $(FFLAGS)  -O 
    CFLAGS := $(CFLAGS)  -O 
  endif
  ifeq ($(MODEL),pio1)
    CPPDEFS := $(CPPDEFS)  -DNO_MPIMOD 
  endif
  ifeq ($(compile_threaded),true)
    LDFLAGS := $(LDFLAGS)  -fopenmp 
  endif
endif
