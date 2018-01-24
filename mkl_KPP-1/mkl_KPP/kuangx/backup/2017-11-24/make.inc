LIBS_FFT = -L /opt/fftw3/lib -lfftw3
LIBS_LAPACK = -llapack -lblas -lpthread
LIBS = $(LIBS_FFT) $(LIBS_LAPACK)
CC := gfortran
MAKEROOT = /home/kuangx/KPP/kuangx/2017-11-24
%.o:%.F90
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS)
