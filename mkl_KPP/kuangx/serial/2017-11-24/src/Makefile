include ../make.inc
OBJ= main.o sparse_run.o sparseGetPlaneSolution.o
SRCS :=$(wildcard *.F90)
OBJS :=$(patsubst %.F90, %.o, $(SRCS))
#sparseGetPlaneSolution.o:sparseGetPlaneSolution.c
#	gfortran -c -o $@ $<
./PHONY:all
all:$(OBJS)
