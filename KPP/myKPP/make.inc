LIBS_FFT = -L /share/home/kuangx/opt/fftw3/lib -lfftw3_mpi -lfftw3
LIBS_LAPACK =-L/share/home/kuangx/opt/lapack/lib -llapack  -lblas  
LIBS_SCALAPACK =-L/share/home/kuangx/opt/scalapack-2.0.2/lib -lscalapack  
LIBS = $(LIBS_FFT) $(LIBS_SCALAPACK) $(LIBS_LAPACK) -lgfortran -lm
CC := mpicc 
MAKEROOT = /home/kuangx/KPP/myKPP
%.o:%.c
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS) 
