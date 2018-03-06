LIBS_FFT = -L/opt/fftw3/lib -lfftw3_mpi -lfftw3
LIBS_SCALAPACK =-L/opt/scalapack/scalapack-2.0.2 -lscalapack -lpthread -lm
LIBS_LAPACK= -llapack -lblas 

LIBS = $(LIBS_FFT) $(LIBS_SCALAPACK) $(LIBS_LAPACK) -lgfortran
CC := mpicc 
MAKEROOT = /home/kuangxiong/Documents/SLEC-4/mkl_KPP-1/mkl_KPP/myKPP
%.o:%.c
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS) 
