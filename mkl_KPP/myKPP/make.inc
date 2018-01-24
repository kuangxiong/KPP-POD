LIBS_FFT = -L /share/home/kuangx/opt/mkl_fftw/lib -lfftw3_mpi -lfftw3
LIBS_SCALAPACK =-L/soft/intel_2018_update1/mkl/lib/intel64_lin/  -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64  -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl 
LIBS = $(LIBS_FFT) $(LIBS_SCALAPACK)  -lgfortran
CC := mpicc 
MAKEROOT = /home/kuangx/mkl_KPP/myKPP
%.o:%.c
	${CC} -c $< -o ${MAKEROOT}/obj/$@ $(LIBS) 
