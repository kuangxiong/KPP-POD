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
//#include"/home/kuangxiong/Software/PHG_source/scalapack-2.0.2/PBLAS/SRC/PBblacs.h"
//#include"/home/kuangxiong/Software/PHG_source/scalapack-2.0.2/PBLAS/SRC/pblas.h"
//#include"/home/kuangxiong/Software/PHG_source/scalapack-2.0.2/SRC/pblas.h"
void
KPP_mpi_svd(int ncup, int N, int localN, complex *data, double *S, complex *U, int nprocs, int myid)
{  
	int context, i, j, k, info, lwork, M, tmpN;
	int descA[9], nprow, npcol,myrow, mycol, mU, nb;
	int  descU[9], mA, ONE=1;
	double  wkopt, *rwork, rwkopt;
    complex *work;

	nprow = nprocs;
	npcol = 1;  
	nb = 1;
	/************************************/
	Cblacs_pinfo(&myid, &nprocs);
	Cblacs_get( -1, 0, &context);
	Cblacs_gridinit(&context, "Row", nprow, npcol);
	Cblacs_gridinfo(context, &nprow, &npcol, &myrow, &mycol);
	/**************************************************************/
//	if(myid==1)
//	for(k=0; k< ncup; k++){
//		for(i=0; i< localN; i++)
//			for(j=0; j< N; j++)
//				printf("%e&%e\t", data[k*(localN*N)+i*N+j]);
//		printf("\n");
//	}
    M = N * N;	
	mA = localN * N;
	mU = localN * N;
	/*******************************************************************/
    descA[0] = 1; descA[1] = 0; descA[2] =M; descA[3]= ncup; descA[4]= 1; 
    descA[5]= 1; descA[6] = 0; descA[7] = 0; descA[8]= mA;
    descU[0] = 1; descU[1] = 0; descU[2] =M; descU[3]= M; descU[4]= 1; 
    descU[5]= 1; descU[6] = 0; descU[7] = 0; descU[8]= mU;
	/********************************************************************* */
	lwork= -1;
	info=0;
	pzgesvd_("V", "N", &M, &ncup, data, &ONE, &ONE, descA, S, U, &ONE, &ONE, descU, NULL,  NULL, NULL,  NULL, &wkopt, &lwork , &rwkopt, &info);

	ONE = 1;
	lwork=(int)wkopt;
    tmpN =(int)rwkopt;
//	printf("hello:%d\t%d\n", lwork, tmpN);
	work = (complex*)malloc(lwork * sizeof(complex));
	rwork = (double*)malloc(tmpN * sizeof(double));
	pzgesvd_("V", "N", &M, &ncup, data, &ONE, &ONE, descA, S, U, &ONE, &ONE, descU, NULL,  NULL, NULL,  NULL, work, &lwork, rwork,  &info);

//	printf("hello\n");
//	Cblacs_gridexit(context);
	Cblacs_gridexit(0);
	free(work); 
	free(rwork);
}
