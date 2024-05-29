#-*-Makefile-*-

COMP = gnu

#FC=gfortran
#FFLAGS=-fdefault-real-8 -fdefault-double-8 -cpp -ffree-line-length-none -fopenmp
#FWARN= #-Wall -fbounds-check
#FOPT=-march=native -flto -funroll-loops -O3 #-fblas  (check this blas one to replace matmul calls)
#FPROF= #-g -pg

ifeq ($(COMP),intel)
	#Intel Compiler
	FC=ifort
	FOPT= -O3 -qopenmp #-xHost -ipo #-parallel
	FFLAGS=-r8 -fpp #-traceback -check all -check noarg_temp_created
	FPROF= #-pg -g
	LIBS = -lfftw3_omp -lfftw3 -mkl -lm
else ifeq ($(COMP),gnu)
	#GNU Compiler
	FC = gfortran
	FOPT = -o2 -fopenmp
	FFLAGS = -fdefault-real-8 -fdefault-double-8 -cpp -ffree-line-length-none -Wall
	FFTW_INC = -I /usr/include
	FFTW_LIB = -L /usr/lib/x86_64-linux-gnu
	BLAS_LIB = -L /usr/lib/x86_64-linux-gnu/blas
	LAPACK_LIB = -L /usr/lib/x86_64-linux-gnu/lapack
	LIBS = -lfftw3_omp -lfftw3 -llapack -lblas -lm
endif

VPATH = 

OBJS = fftw_mod.o params.o vars.o tensor_mod.o grv_mod.o lat_init.o filter_mod.o io.o

spectral: %: $(OBJS) fld_gen_tt.o
	$(FC) $(FFTW_INC) $(FFLAGS) $(FOPT) $(FPROF) -o gen fld_gen_tt.f90 $(OBJS) $(FFTW_LIB) $(BLAS_LIB) $(LAPACK_LIB) $(LIBS)

%.o: %.f90
	$(FC) $(FFTW_INC) $(FFLAGS) $(FOPT) $(FPROF) -c $< -o $@ $(FFTW_LIB) $(BLAS_LIB) $(LAPACK_LIB) $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f gen
