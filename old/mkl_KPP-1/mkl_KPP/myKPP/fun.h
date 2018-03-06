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
#include "/opt/fftw3/include/fftw3-mpi.h"

#define PI 3.1415926535
#define DOUBLE double

static double gamma1 = 0.9999999999;
static double gamma2 = 0.99999;
static double gamma3 = 0.999;
static int intval = 100;  //collect PW solution between two interval time
static double T = 0.52;  //total time
static double tmax = 0.01; // the length for collecting PW solution 
static double errflag = 0.00001; //errindicator threshold
static int initer = 10000;  // number of time for updating POD Matrix

void
zcopy_(int*, complex*, int*, complex*, int*);

void 
KPP_PW(int, DOUBLE, DOUBLE, DOUBLE , DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int , int, MPI_Comm);

void 
testKPP_PW(int, DOUBLE, DOUBLE, DOUBLE , DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int , int, MPI_Comm);

void 
KPP_APOD(int, DOUBLE, DOUBLE, DOUBLE , DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int , int, MPI_Comm);

void 
KPP_TGAPOD(int, DOUBLE, DOUBLE, DOUBLE , DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int , int, MPI_Comm);

void 
KPP_POD(int, DOUBLE, DOUBLE, DOUBLE , DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int , int, MPI_Comm);

void
KPP_ComputePlane(int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, complex*, int, int, long int, long int, MPI_Comm);

void
KPP_BuildPODMatrix(int , int , int , int , int, double , double , 
				   double, double, double, complex *, complex*, complex*, double*, 
				   double *, int , int, MPI_Comm);

void
KPP_BuildPODMatrix1(int , int , int , int , int, double , double , 
				   double, double, double, complex *, complex *, complex *, 
                   complex*, complex*, double*, double *, int , int , MPI_Comm);

double 
KPP_GetErrIndicator(int , int , int , double , double ,complex *, 
                    complex*, complex *, complex *, complex *, MPI_Comm);

void
KPP_Build_FFTMatrix(int, int, int, complex*, int , int, complex*, MPI_Comm);


void
KPP_mpi_svd(int, int, int, complex*, double*, complex *, int, int);

void 
KPP_GetInitialPOD(int, int , int, complex *, complex *, complex *, int, int);

void
KPP_Build_FFTMatrix1(int, int, int, complex*, complex*, MPI_Comm);

int 
KPP_GetPODNumber(double*, int, double);

double 
dnrm2_(int *, complex *, int *);
