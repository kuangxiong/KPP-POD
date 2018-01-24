/*************************************************************************
	> File Name: KPP_SVD_C.c
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2017年12月20日 星期三 22时07分34秒
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"complex.h"
#include"mpi.h"
#include"/home/kuangxiong/Software/PHG_source/scalapack-2.0.2/PBLAS/SRC/PBblacs.h"
//#include"/home/kuangxiong/Software/PHG_source/scalapack-2.0.2/PBLAS/SRC/pblas.h"
//#include"/home/kuangxiong/Software/PHG_source/scalapack-2.0.2/SRC/pblas.h"
void
kpp_mpi_svd_c_(int *ncup1, int *N1, int *localN1, complex data[*ncup1][*localN1*(*N1)], float U[*localN1*(*N1)][*localN1*(*N1)], float S[*ncup1], int *numprocs, int *rank)
{  
    int ncup, N, localN, nprocs, myid;
	int context, i, j, info, lwork, M;
	int descA[9], myrow, mycol, nprow, npcol, mU, nU,nb;
	int  descU[9], mA, nA, ZERO=0, ONE =1;
	double  wkopt, *rwork, rwkopt;
    complex  *work;
    
    ncup = *ncup1;   N = *N1;  localN = *localN1; nprocs=*numprocs;
    myid=*rank; M = N * N; nprow = nprocs;
	npcol = 1;  nb = 1;
	/************************************/
	Cblacs_pinfo(&myid, &nprocs);
	Cblacs_get( -1, 0, &context);
	Cblacs_gridinit(&context, "Row", nprow, npcol);
	Cblacs_gridinfo(context, &nprow, &npcol, &myrow, &mycol);
	/**************************************************************/
	mA = numroc_(&M, &nb, &myrow, &ZERO, &nprow);
	mU = numroc_(&M, &nb, &myrow, &ZERO, &nprow);
	/*******************************************************************/
//	descinit_(descA, &M, &N, &nb, &nb, &ZERO, &ZERO, &context, &mA, &info);
//	descinit_(descU, &M, &M, &nb, &nb, &ZERO, &ZERO, &context, &mU, &info);
    descA[0] = 1; descA[1] = 0; descA[2] =M; descA[3]= ncup; descA[4]= 1; 
    descA[5]= 1; descA[6] = 0; descA[7] = 0; descA[8]= mA;
    descU[0] = 1; descU[1] = 0; descU[2] =M; descU[3]= ncup; descU[4]= 1; 
    descU[5]= 1; descU[6] = 0; descU[7] = 0; descU[8]= mU;
	/********************************************************************* */
	lwork= -1;
	info=0;
	pzgesvd_("V", "N", &M, &ncup, data, &ONE, &ONE, descA, S, U, &ONE, &ONE, &descU, NULL,  NULL, NULL,  NULL, &wkopt, &lwork , &rwkopt, &info);

	lwork=(int)wkopt;
	work = (complex*)malloc(lwork * sizeof(complex));
	rwork = (double*)malloc(rwkopt * sizeof(double));
	
	pzgesvd_("V", "N", &M, &ncup, data, &ONE, &ONE, descA, S, U, &ONE, &ONE, &descU, NULL,  NULL, NULL,  NULL, work, &lwork, rwork,  &info);

	Cblacs_gridexit(0);
	free(work); 
}
