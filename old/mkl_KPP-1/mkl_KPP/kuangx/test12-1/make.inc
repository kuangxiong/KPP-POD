LIBS_FFT = -L/opt/fftw3/lib -lfftw3_mpi -lfftw3
LIBS_LAPACK = -llapack -lblas -lpthread
LIBS = $(LIBS_FFT) $(LIBS_LAPACK)
CC := mpifort 
MAKEROOT =~/KPP/kuangx/test12-1
%.o:%.F90
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS)
