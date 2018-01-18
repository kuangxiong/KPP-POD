LIBS_FFT = -L /opt/fftw/fftw-3/lib -lfftw3
#LIBS_LAPACK = -L /opt/intel/Compiler/11.1/064/mkl/lib/em64t -Wl,-rpath,/opt/intel/Compiler/11.1/064/mkl/lib/em64t -lmkl -lguide -lpthread
#LIBS_LAPACK = -L /opt/intel/mkl/lib/intel64 -Wl,-rpath,/opt/intel/mkl/lib/intel64 -lmkl -lguide -lpthread
LIBS_LAPACK = -llapack -lblas -lpthread

LIBS = $(LIBS_FFT) $(LIBS_LAPACK)
main: main.o sparse_run.o
	@rm -f main
	gfortran -o $@ main.o sparse_run.o  $(LIBS)
#	ifort -o $@ main.o sparse_run.o  $(LIBS)
sparse_run.o: sparse_run.F90
	gfortran -c -o $@ $<
#	ifort -c -o $@ $<
main.o: main.F90
	gfortran -c -o $@ $<
#	ifort -c -o $@ $<
.PHONY:clean
clean:
	rm -rf *.out *sig fort.* *.o
