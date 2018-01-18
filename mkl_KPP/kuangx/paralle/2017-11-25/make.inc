LIBS_FFT = -L /share/home/kuangx/opt/fftw/lib -lfftw3
LIBS_LAPACK = -L/share/home/kuangx/opt/lapack/lib -llapack -lblas -lpthread
LIBS = $(LIBS_FFT) $(LIBS_LAPACK)
CC := gfortran
MAKEROOT = ~/KPP/kuangx/paralle/2017-11-25/
%.o:%.F90
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS)
