LIBS_FFT = -L /share/home/kuangx/opt/fftw/lib -lfftw3
LIBS_LAPACK = -L/share/home/kuangx/opt/lapack/lib -llapack -lblas 
LIBS = $(LIBS_FFT) $(LIBS_LAPACK) -lgfortran -lm
CC := gfortran
MAKEROOT = /home/kuangx/KPP/kuangx/paralle/2017-11-26/
%.o:%.F90
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS)
