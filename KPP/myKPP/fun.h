/*************************************************************************
	> File Name: fun.h
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2017年11月25日 星期六 10时04分51秒
 ************************************************************************/
#ifndef _FUN_H_
#define _FUN_H_
#endif
#include"mpi.h"
#include"math.h"
#include"string.h"
#include"complex.h"
#include"stdlib.h"
#include "/share/home/kuangx/opt/fftw/include/fftw3-mpi.h"

#define PI 3.1415926535
#define DOUBLE double

void
zcopy_(int*, complex*, int*, complex*, int*);

void 
KPP_PW(int, DOUBLE, DOUBLE, DOUBLE , DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int , int, MPI_Comm);

void 
KPP_APOD(int, DOUBLE, DOUBLE, DOUBLE , DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int , int, MPI_Comm);

void
KPP_ComputePlane(int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, complex*, int, int, long int, long int, MPI_Comm);

void
KPP_Build_FFTMatrix(int, int, int, complex*, int , int, complex*);

void
KPP_mpi_svd(int, int, int, complex*, double*, complex *, int, int);

void 
KPP_GetInitialPOD(int, int , int, complex *, complex *, complex *, int, int);

void
KPP_Build_FFTMatrix1(int, int, int, complex*, complex*);

int 
KPP_GetPODNumber(double*, int, double);


