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
KPP_mpi_svd2(int ncup, int N, int localN, complex *data, double *S, complex *U, int nprocs, int myid)
{  
	int context, i, j, k, info, lwork, M, tmpN;
	int descA[9], nprow, npcol,myrow, mycol, mU, nb;
	int  descU[9], mA, ONE=1,descV[9], mV;
	double  wkopt, *rwork, rwkopt;
    complex *work, *V, *tmpdata;

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
    printf("ncup:%d\n", ncup);
    M = N * N;	
	mA = localN * N;
	mU = localN * N;
    mV = ncup/nprow;
    V = calloc(mV*ncup, sizeof(*V));
    tmpdata = calloc(mA*ncup, sizeof(*tmpdata));
    for(i=0; i< mA; i++)
        for(j=0; j< ncup; j++)
            tmpdata[i*ncup+j]=data[i*ncup+j];
	/*******************************************************************/
    descA[0] = 1; descA[1] = 0; descA[2] =M; descA[3]= ncup; descA[4]= 1; 
    descA[5]= 1; descA[6] = 0; descA[7] = 0; descA[8]= mA;
    descU[0] = 1; descU[1] = 0; descU[2] =M; descU[3]= M; descU[4]= 1; 
    descU[5]= 1; descU[6] = 0; descU[7] = 0; descU[8]= mU;
    descV[0] = 1; descV[1] = 0; descV[2] =ncup; descV[3]= ncup; descV[4]= 1; 
    descV[5]= 1; descV[6] = 0; descV[7] = 0; descV[8]= mV;
	/********************************************************************* */
	lwork= -1;
	info=0;
	pzgesvd_("V", "V", &M, &ncup, data, &ONE, &ONE, descA, S, U, &ONE, &ONE, descU, V,  &ONE, &ONE, descV, &wkopt, &lwork , &rwkopt, &info);

	ONE = 1;
	lwork=(int)wkopt;
    tmpN =(int)rwkopt;
//	printf("hello:%d\t%d\n", lwork, tmpN);
	work = (complex*)malloc(lwork * sizeof(complex));
	rwork = (double*)malloc(tmpN * sizeof(double));
	pzgesvd_("V", "V", &M, &ncup, data, &ONE, &ONE, descA, S, U, &ONE, &ONE, descU, V, &ONE, &ONE, descV, work, &lwork, rwork,  &info);

    /***************************************************
     * comfirm  the result
     * **************************************************/
    complex *testU, *testMix, *error, alpha, beta;

    testU = calloc(mA*ncup, sizeof(*testU));

    testMix =calloc(mA*ncup, sizeof(*testMix));

    error =calloc(mA*ncup, sizeof(*error));

    for(i=0; i< mA; i++)
        for(j=0; j<ncup; j++)
            testU[i+j*mA] = U[i+j*mA]*S[j];

    alpha = 1.0;
    beta = 0;
    tmpN = N*N;
    pzgemm_("N", "N", &tmpN,  &ncup, &ncup, &alpha, testU, &ONE, &ONE, descA, V, &ONE, &ONE, descV, 
            &beta, testMix, &ONE, &ONE, descA);

    for(i=0; i< mA; i++)
        for(j=0; j< ncup; j++)
            error[i*ncup+j] = testMix[i*ncup+j]-tmpdata[i*ncup+j];

    if(myid==0)
        for(i=0; i< mA; i++)
        {    for(j=0; j< ncup; j++)
                if(creal(tmpdata[i*ncup+j])>0.00000001||cimag(tmpdata[i*ncup+j])>0.00000001)
                     printf("%e&%e\t%e&%e\n",testMix[i*ncup+j], error[i*ncup+j]);
//            printf("\n");
        }

    /****************************************************/
//	printf("hello\n");
//	Cblacs_gridexit(context);
	Cblacs_gridexit(0);
	free(work); 
	free(rwork);
}
