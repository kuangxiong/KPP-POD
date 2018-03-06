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
KPP_mpi_ev(int ncup, int N, int localN, complex *data, double *S, complex *U, int nprocs, int myid)
{  
	int context, i, j, k, info, lwork, M, tmpN;
	int descA[9], nprow, npcol,myrow, mycol, mU, nb, mA, ONE=1, tmpM, MM;
	double  *a, *work, wkopt , *tmpS;
    complex value, *a_imag, alpha, beta, *tmpU, tmpvalue;

	nprow = nprocs;
	npcol = 1;  
	nb = 1;
    
    tmpM = localN* N;
    MM=N*N;
    tmpN = ncup;
	/************************************/
	Cblacs_pinfo(&myid, &nprocs);
	Cblacs_get( -1, 0, &context);
	Cblacs_gridinit(&context, "Row", nprow, npcol);
	Cblacs_gridinfo(context, &nprow, &npcol, &myrow, &mycol);
	/**************************************************************/
    a = calloc(ncup * ncup, sizeof(*a));
    tmpS = calloc(ncup, sizeof(*tmpS));
    a_imag = calloc(ncup * ncup, sizeof(*a_imag));
    tmpU = calloc(ncup*localN*N, sizeof(*tmpU));
	/*******************************************************************/
    descA[0] = 1; descA[1] = 0; descA[2] =N*N; descA[3]= 1; descA[4]= 1; 
    descA[5]= 1; descA[6] = 0; descA[7] = 0; descA[8]= localN*N;
    /************************************************************************/
/********************************test*******************************/
//    pzdotc_(&MM, &value, &data[1*localN*N], &ONE, &ONE, descA, &ONE, &data[1*localN*N], &ONE, &ONE, descA, &ONE);
//    if(myid==0)
//    printf("value:%f\t%f\n", value);

/*********************************************************************/
    for(i=0; i< ncup; i++){
        for(j=0; j< ncup; j++){
          pzdotc_(&MM, &value, &data[i*localN*N], &ONE, &ONE, descA, &ONE, &data[j*localN*N], &ONE, &ONE, descA, &ONE);
          a[i+j*ncup] = creal(value);
//          if(myid==0)  
//          printf("1111:%f\n,", a[j*ncup+i]);

        }
    }
//    if(myid==1)
//    for(i=0; i< ncup; i++){
//        for(j=0; j< ncup; j++){
//            printf("%f\t", a[j*ncup+i]);
//        }
//        printf("\n");
//    }
//   if(myid==0)
//       for(i=0; i< ncup; i++)
//            printf("%f,", a[i*ncup+i]);
//   if(myid==0)
//       printf("hahahah:%f\n", creal(value));
//    if(myid==1)
//  for(i=1; i< ncup; i++)
//      for(j=0; j< i; j++)
//       printf("matrix:%f \t%f\n", a[i*ncup+j], a[j*ncup+i]);
//#if 1        
	lwork= -1;
	dsyev_("V", "U", &tmpN, a, &tmpN, tmpS, &wkopt , &lwork, &info);

	ONE = 1;
	lwork=(int)wkopt;
//  printf("lwork:%d\n", lwork);
	work = (double *)malloc(lwork * sizeof(double ));
	dsyev_("V", "U", &tmpN, a, &tmpN, tmpS,  work, &lwork, &info);
    
    for(i=0; i< ncup; i++){
        if(tmpS[ncup-1-i]<0)
            tmpS[ncup-1-i]=0;
        S[i] = sqrt(tmpS[ncup-1-i]);
    }
    for(i=0; i< ncup*ncup; i++){
        a_imag[i] = a[i]+0*I;
    }
    MPI_Barrier(MPI_COMM_WORLD);
//    if(myid==1)
//    for(i=0; i< ncup; i++){
//        for(j=0; j< ncup; j++){
//            printf("%f\t", a[j*ncup+i]);
//        }
//        printf("\n");
//    }
//   if(myid==0)
//   for(i=0; i< ncup; i++){
//       for(j=0; j< ncup; j++){
//           printf("%f,%f, %f, %f\n", a[i], a[ncup+i], a[2*ncup+i], a[3*ncup+i]);
//       }
//       printf("\n");
//   }
    alpha = 1.0;
    beta = 0;
    zgemm_("N", "N", &tmpM, &tmpN, &tmpN, &alpha, data, &tmpM , a_imag, &tmpN,  &beta, 
            tmpU, &tmpM);
    for(i=0; i< ncup; i++){
        if(tmpS[i]>0.000001){
            value = 1.0/sqrt(tmpS[i])+0*I;   
            zaxpy_(&tmpM, &value, &tmpU[i*tmpM], &ONE, &U[i*tmpM], &ONE);
        }
    }
//    pzdotc_(&MM, &value, &U[2*localN*N], &ONE, &ONE, descA, &ONE, &U[2*localN*N], &ONE, &ONE, descA, &ONE);
//    printf("value:%f&%f\n", value);
//    if(myid==0)
//        for(i=0; i< localN*N; i++)
//            printf("U:%f&%f\t%f&%f\n",data[i], U[i]);

//    Cblacs_gridexit(context);
    Cblacs_gridexit(0);
	free(work);
//#endif
}
