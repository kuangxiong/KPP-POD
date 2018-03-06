LIBS_FFT = -L /opt/fftw/fftw-3/lib -lfftw3
LIBS_LAPACK = -llapack -lblas -lpthread
LIBS = $(LIBS_FFT) $(LIBS_LAPACK)
CC := gfortran
MAKEROOT = /home/kuangx/KPP/kuangx/serial/2017-11-23
%.o:%.F90
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS)
