/*************************************************************************
	> File Name: KPP_SVD_C.c
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2017年12月19日 星期二 20时33分26秒
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"mpi.h"
#include"complex.h" 


void kpp_build_fftmatrix_c_(int *n1, int *N1, int *localN1, complex mpiUP_hat[*n1][*localN1*(*N1/2+1)], 
                    int *myrank1, int *nprocs1, complex tmp_hat1[*n1][*localN1*(*N1)])
{
    int i, j, k, dist1, dist2;
    int ncup, N, localN, myrank, nprocs;
    complex tmp_v[*n1][*localN1 * (*N1/2-1)];
    MPI_Status status;

    ncup =*n1; N=*N1; localN = *localN1; myrank=*myrank1; nprocs=*nprocs1;

    //build matrix tmp_hat1  for SVD decomposition
    for(k = 0; k < ncup; k++){
        for(i=0; i< localN; i++){
            for(j=0; j<= N/2; j++){
                tmp_hat1[k][j + i * N] = mpiUP_hat[k][j + i * (N/2+1)];   
            }
        }
    }
    //build matrix tmp_v prepare for tmp_hat1 
    for(k=0; k< ncup; k++){
        for(i=0; i< localN; i++){
            for(j=1; j< N/2; j++){
                tmp_v[k][j + i * (N/2-1)] = conj(mpiUP_hat[k][N/2 - j + i * (N/2+1)]);   
            }
        }
    }
    // get target process for echo rank 
    dist1 = nprocs - myrank -1;
    dist2 = nprocs - myrank;
    if(myrank==0||myrank==nprocs/2)
        dist2 = myrank; 

//    printf("testrank:%d\t%d\n", myrank, dist2);
    
    for(k=0; k<ncup; k++){
        for(i = 1; i<localN; i++){
          MPI_Sendrecv(&tmp_v[k][0+i * (N/2-1)], N/2-1, MPI_COMPLEX, dist1, 991,
                       &tmp_hat1[k][N/2 + 1 + (localN - i) * N], N/2-1, MPI_COMPLEX,
                       dist1, 991, MPI_COMM_WORLD, &status);
        }
    }
    
    for(k=0; k<ncup; k++){
         MPI_Sendrecv(&tmp_v[k][0], N/2-1, MPI_COMPLEX, dist2, 992,
                      &tmp_hat1[k][0], N/2-1, MPI_COMPLEX,
                     dist2, 992, MPI_COMM_WORLD, &status);
    }

    if(myrank==0 || myrank==nprocs/2){
        for(k=0; k< ncup; k++){
            for(i=1; i< N/2; i++){
                tmp_hat1[k][N/2+i] = conj(tmp_hat1[k][N/2-i]);
            }
        }
    }
}

