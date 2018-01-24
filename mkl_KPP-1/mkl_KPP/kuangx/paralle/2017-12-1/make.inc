LIBS_FFT = -L /opt/fftw/fftw-3/lib -lfftw3_mpi -lfftw3 
LIBS_LAPACK = -L/opt/scalapack/lib/ -lscalapack -llapack -lblas -lpthread
LIBS = $(LIBS_FFT) $(LIBS_LAPACK)
CC := mpifort 
MAKEROOT =~/KPP/kuangx/paralle/2017-12-1
%.o:%.F90
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS)
