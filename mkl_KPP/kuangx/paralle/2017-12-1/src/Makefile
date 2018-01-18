include ../make.inc
SRCS :=$(wildcard *.F90)
OBJS :=$(patsubst %.F90, %.o, $(SRCS))
OBJS1 :=KPP_SVD_C.o KPP_Build_FFTMatrix_C.o
#sparseGetPlaneSolution.o:sparseGetPlaneSolution.c
#	gfortran -c -o $@ $<
KPP_Build_FFTMatrix_C.o:KPP_Build_FFTMatrix_C.c
	mpicc -c KPP_Build_FFTMatrix_C.c -o ../obj/KPP_Build_FFTMatrix_C.o
KPP_SVD_C.o:KPP_SVD_C.c
	mpicc -c KPP_SVD_C.c -o ../obj/KPP_SVD_C.o
./PHONY:all
all:$(OBJS) $(OBJS1)
